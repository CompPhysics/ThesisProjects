# Paper tasks — working checklist

Companion to `article-plan.md`. Mark items `[x]` as done. Re-order freely. Notes after each item say what it unblocks or depends on.

---

## 0. Pre-decisions (do these first, in one sitting)

- [ ] **P1. ("(F)" in title) — DEFERRED.** Decide after Phase A + B land, based on remaining time budget for Phase C. Phase C is now explicitly conditional.
- [x] **P2. FCI reference code available.** PRA-era code works. Unblocks A4, A10, C2, C3. *Use split by regime (V1 finding): static/ITE refs (A4, C2, C3) use it as-is; the time-dependent ref (A10) needs the propagation timestep converged first. Heidelberg cross-check dropped.*
- [x] **P3. Decoherence paragraph: IN.** One paragraph in Sec. 7 + Appendix F proof-of-concept. See Phase D below (sits between Phase C and the writing pass).

---

## 1. De-risk the stack (Phase A, validation first)

- [x] **V1. TD-FCI cross-check on the 2-qubit reduction (A10).** Configure the 2-qubit reduction of the current 4-well code; compare the MCTDH(F) **√iSWAP** fidelity against *converged* TD-FCI on the same system. *Load-bearing: if this disagrees, everything downstream is suspect.* Depends on P2. Result recorded in `article-plan.md` (A10).
- [ ] **V2. Drift quick-look from saved `y(t)`.** Notebook cell or short helper that loads saved checkpoints and plots E(t), ‖ψ‖², (⟨S²⟩ where applicable) against the §2.6 tolerances. Likely partly there in `rte_step_plot.ipynb`.
- [x] ~~V3. Wall-time + peak-memory instrumentation.~~ *Dropped — scaling argument now analytical (§2.8); see A9.*

---

## 2. Phase A — motional qubits (mostly executing existing plan)

- [x] **A1. Stage-1 ITE anchors verified converged.** Resolution + reframed anchor-quality criterion in `article-plan.md` (§4.2, §2.6).
- [x] **A2. Stages 2/3/4 anchors converged.** Partner fidelities balanced ~0.54 ≥ 0.49 (verification pass 2026-06-03). Canonical: `chain_step1.npz` + `chain_global_20260512_181719_step{2,3,4}_*`.
- [x] **A3. Stage 5 (long-range w0w3; now Gate-2b direct-arm evidence) feasibility — FEASIBLE.** Re-optimized with `alpha_ordering = "none"`; balanced w0w3 superposition, partner fids 0.498/0.493 (min 0.4929 ≥ 0.49), partner eigenstates near-degenerate (gap 1.7e-3). Canonical: `ITE_trajectories/reopt_step5_from1_stage5_superpos_w0w3.npz`.
- [x] **A4. Convergence study, 2-qubit reduction.** SPF ∈ {2,3,4,5} vs FCI, at the correct optimized geometry. DVR N=64 converged (machine precision). Static energy converges 2→3→4 (9e-2→5e-4→1e-6), plateau at 5. **Dynamic gate (x64): per-state 3 SPF = 0.99999 (0.991 complex64), 4 SPF = 0.99998, 5 SPF over-completes → 0.75; SA needs 5 SPF (0.75 at 3,4 → 1.0 at 5).** Production per-state = 3 SPF; **SA-mode would need 5 (SA@4 insufficient)**. eps=1e-6 finite throughout (the earlier "4+ NaN / fragile" was a wrong-geometry artifact, retracted). Writeup `A4_DVR_convergence/A4_results.md`.
- [x] **A4b. SPF transferability check at production SPF.** 4-qubit anchor DONE (`run_4well_spf.py`, config-1 idle self-convergence): weighted states converged by norb=3–4, norb=4 well-conditioned (min occ 3e-5). Gate-dynamics half CLOSED (A5 Runs 2b–2d): per-state norb=4 clean through the full leg at dt=5e-3, all §2.6 gates pass, spectators exonerated per-state. 3-qubit anchor check dropped (production is 4-well; see Appendix A wording).
- [ ] **A5. Leg-primitive BOBYQA.** *(Gate restructure 2026-06-10: Protocols A/B/C dropped; Section 4 = Gate 1 headline adjacent √iSWAP + Gate 2 same-gate comparison √iSWAP(0,3) routed-vs-direct; full SWAP survives only as routing primitive.)* Optimize per-leg, 2-D over (τ_ramp, τ_hold), 4-state condensed loss with switchable M_ref: (i) adjacent √iSWAP leg (warm start τ_hold ≈ 12–13); (ii) two full-SWAP routing legs (~2× hold). Save traces; verify improvement over adiabatic guess. *FoM PINNED (§2.6):* basis = ITE anchors by dominant configuration, never energy ordering. *Status 2026-06-11 (post Runs 3/4):* both legs CONVERGED **on the wrong loss** (1 − F_xc, leakage-blind — §2.6 metric corrected same day): √iSWAP F_xc 0.926 @ (44.4, 14.1), SWAP F_xc 0.883 @ (15.0, 19.5); leakage-penalized F̄ ≈ 0.78/0.74, **|11⟩ ~64% out of subspace at norb=4 on both**. Optima kept as warm starts; τ_ramp [15,45] clipped at opposite ends (√iSWAP upper, SWAP lower) → widen before re-opt. Log: `A5_gate_leg/A5_run_log.md`.
- [ ] **A5c. Solver-caps spot-check.** One |11⟩ propagation at dt=5e-3 with the original caps (min_tol 1e-5, max_iter 20); assert ‖ψ‖² drift ≤ 1e-6 and time it vs the existing 1e-8/40 run (~12 min). Gates reverting the diagnostic-era caps for subsequent launches (A5d/A6+); does NOT gate Runs 3/4 (loose ≡ tight in M, A5 Run 2c).
- [x] **A5d. Re-optimise both legs with the corrected loss 1 − F̄₄.** Superseded by the global sweep_map (below): BOBYQA polishing moot once the landscape was mapped.
- [~] **A5-sweep. Global (τ_ramp, τ_hold) F̄ map (`sweep_map.py`) — DONE; strategic decision PENDING (user, 2026-06-16 AM).** Two full maps (τ_ramp 2–70) prove the √iSWAP F̄ ceiling is **~0.85, control-limited** — the |11>→{|20>,|02>} doubles excursion caps it at every ramp speed (short/diabatic ramps WORSE, 0.76–0.82; closes only adiabatically = identity). **norb-5 refuted** (basis captures the excursion, sum~0.99). Swap routing leg revives at short ramps (F̄~0.77). Full result + diagnostics: `A5_gate_leg/A5_run_log.md` (Runs sweep + short-ramp); memory `project_motional_sqiswap_ceiling`. *Decision needed:* (a) richer control DOF to close |11>; (b) accept ~0.85 as a reported physical ceiling; (c) lean on §5 spin CPHASE for the high-fidelity entangler.
- [ ] **A6. Direct long-range √iSWAP(0,3) leg (Gate 2b).** BOBYQA at the long-range geometry; A3 anchor = feasibility evidence; gate time is itself a result (report honestly if the hold exceeds the feasible window).
- [ ] **A7. Compose + report K=16.** Gate 1 vs √iSWAP⊗I; Gate 2a routed 5-leg composition and Gate 2b direct, both vs the same ideal √iSWAP(0,3)⊗I. Report F̄₁₆ per arm (F_xc shape diagnostic only), leakage, total-T ratio, per-leg error accumulation.
- [ ] **A8. Sensitivity sweep.** ±5% global scaling of x(t) on optimised Gate 1 + Gate 2 arms.
- [ ] **A8b. eps spot-check.** Re-run one optimised leg (Gate-1 √iSWAP) at eps 1e-5 vs 1e-6; report ΔF in the §2.6 error budget. (Mixed-precision half dropped — production RTE is full x64 on CPU per §2.4 / simulator skill.) One leg, hours.
- [x] **A9. Analytical scaling crossover numbers.** All seven §2.8 numbers confirmed (within 11% of exact; fine at 1 s.f.). Derivations + table in `A9_scaling_crossover/A9_results.md` (script `check_numbers.py`). Conventions: complex128 16 B/amp, decimal units, state-vector-only (conservative). *Wording flag resolved:* §2.8 "grows polynomially" replaced with the precise statement (orbitals linear, A-tensor small-base 4^K vs 64^K; polynomial route for large K = ML-MCTDH layering).
- [ ] **A11. Phase-A writing pass.** Section 4 of paper, including the figures it owns. *(Doing this before Phase B forces honesty about what Phase A actually shows. A10 = TD-FCI cross-check, tracked as V1 above — number kept free to match article-plan.)*

---

## 3. Phase B — spin qubits

- [ ] **B0. Decide go/no-go on Phase B.** Based on whether Phase A landed cleanly + time budget. If no-go: drop Sec. 5, retitle, re-scope. *(Honest checkpoint, not a formality.)*
- [ ] **B0.5. Analytic two-qubit coupling estimate.** Differential Coulomb shift between double-dot charge configurations → implied CPHASE gate time. Pencil-and-paper, gates B6/B7: if gate time ≫ feasible propagation window, invoke the §5.7 fallback before building the 4-species pipeline.
- [ ] **B1. Build double-dot geometry module.** Global DVR per double-dot, parallel to `optimize_qubit_params.py`.
- [ ] **B2. Spin-as-species ITE.** Two species (↑, ↓), 1 e each. Loss modules for |S⟩, |T₀⟩ anchors at idle + detuned configurations, paralleling `losses/stage*.py`.
- [ ] **B3. (1,1)–(0,2) SPF convergence sweep.** At most strongly-detuned anchor. Charge populations, ⟨S²⟩, compare to Hubbard-limit reference. Pin production SPF.
- [ ] **B4. ε(t)-schedule single-qubit forward propagation.** Project onto |S⟩, |T₀⟩ at saved t. Verify ⟨S²⟩ drift < 10⁻³.
- [ ] **B5. Single-qubit gate BOBYQA.** Save trace, report improvement.
- [ ] **B6. Two-qubit geometry: 4 species.** ↑_A, ↓_A, ↑_B, ↓_B on full 4-well.
- [ ] **B7. Two-qubit gate BOBYQA.** *Slip risk; fallback per §5.7 is to drop B6/B7 and present single-qubit only.*
- [ ] **B8. Phase-B writing pass.** Section 5.

---

## 4. Phase C — multi-electron-per-dot (earns the (F))

- [ ] **C0. Decide go/no-go on Phase C.** Tied to P1.
- [ ] **C1. Intra-species `g_dict` via `SineDVR.get_g`.** Soft-Coulomb with ε > 0 regulariser for r₁ = r₂. Unit-style check that `SemiDirectCI` N_s = 2 path works.
- [ ] **C2. Isolated 2-e well ITE.** Compare ground-state energy + NO occupations to FCI on that well. Target: ≤ 10⁻³ energy gap. Depends on P2.
- [ ] **C3. Coupled 2e + 1e benchmark.** vs small-basis 3-electron FCI on 2-well subsystem. *This is the "(F)" earner.*
- [ ] **C4. Short forward propagation** with one multi-electron well under a Protocol A leg. One figure.
- [ ] **C5. Phase-C writing pass.** Section 6.

---

## 5. Phase D — decoherence proof-of-concept (Appendix F)

Single-oscillator bath coupled to one motional qubit. Small scope by design; honest about revivals and zero-T.

- [ ] **D1. Bath species setup.** One species, 1 particle, parabolic potential on a sine-/HO-DVR. Bath frequency ω_B chosen so 2π/ω_B exceeds the gate window.
- [ ] **D2. Bilinear coupling.** `g_dict[(qubit, bath)] = κ · x_qubit · x_bath`. Pick κ to give a visible coherence drop on the gate timescale.
- [ ] **D3. Bath initial state.** Ground state via ITE; coherent-state seed as sanity check.
- [ ] **D4. Reduced density matrix routine.** Trace out the bath species at saved t. Reuse 1-RDM machinery in `SemiDirectCI_n_species` where possible; otherwise a short helper.
- [ ] **D5. Run + plot.** Propagate one motional qubit + bath through an idle window and one Protocol-A leg. Plot ρ_S off-diagonal decay vs κ.
- [ ] **D6. Sec. 7 decoherence paragraph.** Frames the bath-as-species idea, points to Appendix F, cites Wang/Thoss and standard system-bath MCTDH literature.
- [ ] **D7. Appendix F.** Single figure, ½–1 page of text. Frame as entanglement-induced coherence loss into a single monitored mode; honest about revivals + zero-T limitation.

---

## 6. Writing & wrap-up

- [ ] **W1. Methods section (Sec. 2)** with pinned fidelity formula, error budget, scope-disclaimer paragraph.
- [ ] **W2. Model section (Sec. 3).** Generic 1-D quantum-dot model, platform-agnostic framing.
- [ ] **W3. Discussion (Sec. 7).** One paragraph each on: encoding generality, scaling, limits, comparison to MCTDH-X / ML-MCTDH(F) / TD-DMRG.
- [ ] **W4. Conclusion (Sec. 8).**
- [ ] **W5. Appendix A** (convergence + transferability).
- [ ] **W6. Appendix B** (spin-as-species formalism).
- [ ] **W7. Appendix C** (numerical details).
- [ ] **W8. Appendix D** (failed autodiff QOCT — short, honest).
- [ ] **W9. Appendix E** (reproducibility).
- [ ] **W10. Final figure pass.** ≤ 7 numbered figures (matches the article-plan figure list), captions written.
- [ ] **W11. End-to-end read-through** with the verification criteria from `article-plan.md` open alongside.

---

## Standing reminders

- Production RTE = full x64 VMF on CPU, **per-state, norb=4, dt=5e-3** (`.claude/skills/simulator/SKILL.md` is the run-settings source of truth; SA infeasible with idle spectators — §2.4). Complex64/x32 is the GPU/CMF fast path only; if used: always `orbital_x64=True` + `e_shift="auto"` + builder downcasts. (See memory.)
- `rte_step.py` default DT=1e-2 sits at the Newton stability edge (~7% transfer error, norm drift to 2e-3 on correlated states). Don't reuse dt=1e-2-era numbers as converged; production dt=5e-3.
- Audit diffs before any multi-hour run.
- Re-evaluate this list after each phase boundary, not just at the end.
