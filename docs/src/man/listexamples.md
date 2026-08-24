# List of examples
Examples with the letter A refer to diffusion models in a single material, while examples with the letters B and C refer to a diffusion couple. In B the model uses flux conditions to describe the interface whereas the models labelled C perform total mass balance calculations. Examples with letter D describe the thermodynamically constrained crystal growth/resorption. However, we have some more examples included, which show the stepwise increasing complexity of the codes. They have names without any reference letter.

## Examples

1. A1: Intracrystalline diffusion within a planar crystal
2. A2: Intracrystalline diffusion in a spherical crystal
3. B1: Intercrystalline diffusion within a spherical diffusion couple
4. B2: Diffusion within a spherical diffusion couple for the case of time-evolving  diffusivity
5. B3: Major element diffusion within a diffusion couple    => **Pending**
6. B4: Spherical crystal growth due to Rayleigh fractionation in a growth and diffusion couple with ``D^A \ll D^B``
7. B5: Growth of an alloy from a melt in a planar geometry
8. B6: Growth of a spherical crystal in a diffusion couple (``v_A > 0``)
9. B7: Growth of a spherical crystal in a diffusion couple (``v_A < 0``)
10. C1: Spherical crystal growth due to Rayleigh fractionation in a growth and diffusion couple with ``D^A \ll D^B``
11. C2: Growth of a spherical crystal in a diffusion couple (``v_A < 0``)
12. D1: Diffusion-limited crystal growth of olivine
13. D2: Diffusion-limited crystal growth of olivine with a user-defined, non-monotonic  temperature-time path   => **Pending**
14. `Simple_Diff`: Diffusion within a single crystal
15. `Diff_couple_no_interaction`: Diffusion couple, which is built from 2 single crystals without an ion-exchange reaction
16. `Diff_couple_Flux`: Diffusion couple with ion-exchange using flux balance at the interface
17. `Diff_couple_MB`: Diffusion couple with ion-exchange using total mass balance at the interface
18. `Diff_couple_Flux_growth`: Diffusion couple with ion-exchange using flux balance at the interface and simultaneous growth
19. `Diff_couple_MB_growth`: Diffusion couple with ion-exchange using total mass balance at the interface and simultaneous growth

## Correspondence to the companion paper

The benchmarked examples (A, B, C, D) reproduce the figures of [Stroh2025](@cite) one-to-one:

| Example | Figure in [Stroh2025](@cite) |
|:-------:|:------|
| A1      | Fig. 2 |
| A2      | Fig. 3 |
| B1      | Fig. 4 |
| B2      | Fig. 5 |
| B4      | Fig. 6 |
| B5      | Fig. 7 |
| B6      | Fig. 8 |
| B7      | Fig. 9 |
| C1      | Fig. C1 (Appendix C; same setup as B4/Fig. 6, but with the total mass-balance condition) |
| C2      | Fig. C2 (Appendix C; same setup as B7/Fig. 9, but with the total mass-balance condition) |
| D1      | Fig. 10 |

B3 has no figure (it is not yet implemented, see above) and D2 is a documentation example added after publication, not part of the paper.
