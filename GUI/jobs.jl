using UUIDs

"""Thrown by a backend's time loop when `should_stop()` returns true, so a cancelled run is
distinguishable from a genuine error."""
struct SimulationCancelled <: Exception end

"""
A single simulation run, tracked from `:running` through `:done`/`:failed`/`:cancelled`.
`result_dir` holds the name of the folder under `GUI/results/` containing the saved
plot/data/inputs once done; `summary` is a one-line human-readable result summary;
`error` holds a readable message if the run failed; `cancel_requested` is set by the
`/jobs/{id}/cancel` endpoint and checked by the running backend function itself.
"""
mutable struct Job
    id::String
    mode::Symbol
    task::Task
    status::Symbol            # :running, :done, :failed, :cancelled
    result_dir::Union{String,Nothing}
    summary::Union{String,Nothing}
    error::Union{String,Nothing}
    cancel_requested::Bool
    started_at::Float64        # time(), used for elapsed-time display
end

# This GUI is a local, single-user tool: only one simulation runs at a time,
# tracked in a single global slot rather than a general job queue.
const CURRENT_JOB = Ref{Union{Job,Nothing}}(nothing)
const JOB_LOCK = ReentrantLock()

"""
    start_job!(f, mode) -> job_id::String or nothing

Starts `f` (a one-argument closure `f(should_stop)` that returns `(result_dir, summary_string)`,
see `save_run_outputs` in `results.jl`) on a background thread via `Threads.@spawn`, so it
doesn't block the HTTP server thread. `f` is called with a zero-argument `should_stop`
function it should pass through to the backend's `should_stop` keyword argument. Returns the
new job's id, or `nothing` if a job is already running (caller should respond with a "please
wait" message in that case). `f` is the first argument (rather than `mode`) so this can be
called with `do...end` block syntax: `start_job!(mode) do should_stop ... end`.
"""
function start_job!(f, mode::Symbol)
    lock(JOB_LOCK) do
        if CURRENT_JOB[] !== nothing && CURRENT_JOB[].status == :running
            return nothing
        end
        id = string(uuid4())
        job = Job(id, mode, Task(() -> nothing), :running, nothing, nothing, nothing, false, time())
        should_stop = () -> job.cancel_requested
        job.task = Threads.@spawn begin
            try
                result_dir, summary = f(should_stop)
                job.result_dir = result_dir
                job.summary = summary
                job.status = :done
            catch e
                if e isa SimulationCancelled
                    job.status = :cancelled
                else
                    job.error = sprint(showerror, e)
                    job.status = :failed
                end
            end
        end
        CURRENT_JOB[] = job
        return id
    end
end

"""Returns the job with the given id, or `nothing` if it doesn't match the current job."""
function get_job(id::AbstractString)
    job = CURRENT_JOB[]
    return (job !== nothing && job.id == id) ? job : nothing
end

"""Requests cancellation of the given job, if it's the current, still-running job. Returns
`true` if a cancellation was requested, `false` if there was nothing running to cancel."""
function cancel_job!(id::AbstractString)
    job = get_job(id)
    (job === nothing || job.status != :running) && return false
    job.cancel_requested = true
    return true
end
