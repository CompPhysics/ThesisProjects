"""
Finite-N Lipkin model TD dynamics + Parametric Matrix Model (PMM) sequence emulator.

Features:
  1. Finite-N Lipkin Hamiltonian in explicit collective-spin basis (dimension N+1):
        H = epsilon * Jz - (chi / (N - 1)) * Jx^2
     where Jx, Jy, Jz are standard angular momentum operators for total spin J = N/2.

  2. Exact many-body time evolution:
        i d|psi>/dt = H |psi>
     implemented with 4th-order Runge–Kutta for a time-independent Hamiltonian.

  3. Energy conservation diagnostics:
     E(t) = <psi(t)|H|psi(t)>
     Plots E(t) and relative deviation from E(0).

  4. Parametric Matrix Model (PMM) as a sequence model:
     Input at each step:
         x = [sx(t), sy(t), sz(t), epsilon, chi, dt]
     Output:
         y = [sx(t+dt), sy(t+dt), sz(t+dt)]
     The PMM:
         - Uses complex Hermitian primary matrices H(x),
         - Interpreted as generators of unitaries (we only need eigenvectors),
         - Uses complex Hermitian secondary matrices to form outputs.

  5. Emulation of entire trajectories by iterating the learned PMM step:
     Compare exact Bloch components vs PMM-predicted ones.

Dependencies: numpy, torch, matplotlib

Run:
    python lipkin_tdhf_pmm_sequence.py
"""

import math
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt


# ============================================================
# 1. Finite-N Lipkin model in collective-spin basis
# ============================================================

class LipkinQuantum:
    """
    Finite-N Lipkin model in the J-basis (collective spin representation).

    We use J = N/2, dimension d = 2J+1 = N+1, basis |m> with m = -J, ..., J.

    Operators:
        Jz |m> = m |m>
        J+ |m> = sqrt(J(J+1) - m(m+1)) |m+1>
        J- |m> = sqrt(J(J+1) - m(m-1)) |m-1>
        Jx = (J+ + J-) / 2
        Jy = (J+ - J-) / (2i)

    Hamiltonian (time-independent):
        H = epsilon * Jz - (chi / (N - 1)) * Jx^2

    We solve i d|psi>/dt = H |psi>.
    """

    def __init__(self, N=10, epsilon=1.0, chi=1.0):
        assert N >= 2, "N must be >= 2"
        self.N = N
        self.J = N / 2.0
        self.epsilon = epsilon
        self.chi = chi

        self.dim = int(2 * self.J + 1)
        self._build_operators()
        self._build_hamiltonian()

    def _build_operators(self):
        J = self.J
        dim = self.dim

        # Basis index: m = -J,...,J mapped to 0,...,dim-1
        ms = np.arange(-J, J + 1, 1, dtype=float)

        Jz = np.diag(ms)
        Jp = np.zeros((dim, dim), dtype=np.complex128)
        Jm = np.zeros((dim, dim), dtype=np.complex128)

        for i, m in enumerate(ms):
            # Raising: m -> m+1
            if i < dim - 1:
                mp = m + 1
                coeff = math.sqrt(J * (J + 1.0) - m * mp)
                Jp[i + 1, i] = coeff

            # Lowering: m -> m-1
            if i > 0:
                mm = m - 1
                coeff = math.sqrt(J * (J + 1.0) - m * mm)
                Jm[i - 1, i] = coeff

        Jx = 0.5 * (Jp + Jm)
        Jy = (Jp - Jm) / (2.0j)

        self.Jz = Jz
        self.Jx = Jx
        self.Jy = Jy
        self.Jp = Jp
        self.Jm = Jm

    def _build_hamiltonian(self):
        N = self.N
        eps = self.epsilon
        chi = self.chi

        Jz = self.Jz
        Jx = self.Jx

        # Lipkin Hamiltonian:
        #   H = eps * Jz - (chi / (N - 1)) * Jx^2
        # N-1 in denominator is one common convention; you can adjust as needed.
        H = eps * Jz - (chi / (N - 1.0)) * (Jx @ Jx)
        # Ensure Hermiticity numerically
        H = 0.5 * (H + H.conj().T)

        self.H = H

    def coherent_state(self, theta, phi):
        """
        Construct a spin-coherent state |theta, phi> in the J-basis.

        |theta,phi> = sum_m C_m |m>, with
            C_m = sqrt(C(2J, J+m)) (cos(theta/2))^{J+m}
                                      (sin(theta/2) e^{i phi})^{J-m}
        """
        J = self.J
        dim = self.dim
        ms = np.arange(-J, J + 1, 1, dtype=float)

        coefficients = np.zeros(dim, dtype=np.complex128)
        for idx, m in enumerate(ms):
            k = int(J + m)  # 0,...,2J
            n = int(2 * J)
            # binomial coefficient C(n, k)
            C = math.comb(n, k)
            coeff = (math.cos(theta / 2.0) ** (J + m) *
                     (math.sin(theta / 2.0) ** (J - m)) *
                     np.exp(1j * (J - m) * phi))
            coefficients[idx] = math.sqrt(C) * coeff

        # Normalize
        norm = np.linalg.norm(coefficients)
        if norm > 0:
            coefficients /= norm

        return coefficients

    def rhs(self, psi):
        """
        Right-hand side of Schrödinger equation:
            i d|psi>/dt = H |psi|
        => d|psi>/dt = -i H |psi|
        """
        return -1.0j * self.H @ psi

    def evolve_rk4(self, psi0, t_max, n_steps):
        """
        Evolve |psi> using 4th-order Runge–Kutta.

        Returns:
            times: (n_steps+1,)
            psis : (n_steps+1, dim)
        """
        dt = t_max / n_steps
        times = np.linspace(0.0, t_max, n_steps + 1)
        psis = np.zeros((n_steps + 1, self.dim), dtype=np.complex128)
        psis[0] = psi0

        psi = psi0.copy()
        for n in range(n_steps):
            k1 = self.rhs(psi)
            k2 = self.rhs(psi + 0.5 * dt * k1)
            k3 = self.rhs(psi + 0.5 * dt * k2)
            k4 = self.rhs(psi + dt * k3)

            psi = psi + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

            # Normalize to avoid numerical drift
            norm = np.linalg.norm(psi)
            if norm > 1e-14:
                psi = psi / norm

            psis[n+1] = psi

        return times, psis

    def expectation_J(self, psi):
        """
        Compute Bloch-like components:
            sx = <Jx>/J, sy = <Jy>/J, sz = <Jz>/J
        for state |psi>.
        """
        psi_dag = psi.conj().T
        Jx = self.Jx
        Jy = self.Jy
        Jz = self.Jz
        J = self.J

        sx = np.real(psi_dag @ (Jx @ psi)) / J
        sy = np.real(psi_dag @ (Jy @ psi)) / J
        sz = np.real(psi_dag @ (Jz @ psi)) / J
        return np.array([sx, sy, sz], dtype=np.float64)

    def expectation_energy(self, psi):
        """
        Compute E = <psi|H|psi>.
        """
        psi_dag = psi.conj().T
        return np.real(psi_dag @ (self.H @ psi))


# ============================================================
# 2. Dataset generation for PMM sequence learning
# ============================================================

def generate_sequence_dataset(num_trajectories=50,
                              N=10,
                              n_steps=200,
                              t_max=10.0,
                              epsilon=1.0,
                              chi_min=0.1,
                              chi_max=1.5,
                              subsample=1,
                              seed=1234):
    """
    Generate a dataset of one-step transitions from exact Lipkin dynamics.

    For each trajectory:
      - Sample interaction chi in [chi_min, chi_max].
      - Sample initial coherent state angles (theta0, phi0).
      - Evolve quantum state with LipkinQuantum.
      - Extract Bloch-like vector s(t) = <J>/J.
      - For each step n -> n+1 (optionally subsampled) form:
            input x_n  = [s_x(t_n), s_y(t_n), s_z(t_n), epsilon, chi, dt/t_max]
            target y_n = [s_x(t_{n+1}), s_y(t_{n+1}), s_z(t_{n+1})]

    Returns:
        X: numpy array (N_pairs, 6)
        Y: numpy array (N_pairs, 3)
    """
    rng = np.random.default_rng(seed)
    dt = t_max / n_steps

    X_list = []
    Y_list = []

    for traj in range(num_trajectories):
        chi = rng.uniform(chi_min, chi_max)

        # Random initial coherent state
        theta0 = np.arccos(2 * rng.random() - 1.0)  # uniform on sphere
        phi0 = 2 * np.pi * rng.random()

        model = LipkinQuantum(N=N, epsilon=epsilon, chi=chi)
        psi0 = model.coherent_state(theta0, phi0)
        times, psis = model.evolve_rk4(psi0, t_max=t_max, n_steps=n_steps)

        # Precompute Bloch vectors
        S = np.zeros((n_steps + 1, 3), dtype=np.float64)
        for i in range(n_steps + 1):
            S[i] = model.expectation_J(psis[i])

        for n in range(0, n_steps, subsample):
            s_n = S[n]
            s_np1 = S[n + 1]  # we must have n_steps >= n+1

            x_n = np.array([
                s_n[0],
                s_n[1],
                s_n[2],
                epsilon,
                chi,
                dt / t_max
            ], dtype=np.float64)

            X_list.append(x_n)
            Y_list.append(s_np1)

    X = np.stack(X_list, axis=0)
    Y = np.stack(Y_list, axis=0)
    return X, Y


# ============================================================
# 3. Unitary-style PMM (complex Hermitian primary + secondary)
# ============================================================

class UnitaryPMMSequence(nn.Module):
    """
    PMM for one-step map:
        x -> s_next

    x has dimension p, s_next has dimension q.
    Primary matrices: complex Hermitian
        H(x) = H0 + sum_l x_l H_l
    Secondary matrices: complex Hermitian S_k.

    We use eigenvectors of H(x) and form:
        y_k(x) = bias_k + sum_{i=1}^r Re( v_i^\dagger S_k v_i )

    where v_i are eigenvectors associated with the lowest r eigenvalues.
    """

    def __init__(self, input_dim, output_dim, latent_dim=8, r=3):
        super().__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.latent_dim = latent_dim
        self.r = r

        # Complex parameters; will be symmetrized to Hermitian.
        self.H0 = nn.Parameter(
            0.1 * (torch.randn(latent_dim, latent_dim, dtype=torch.complex64))
        )
        self.H_inputs = nn.Parameter(
            0.1 * (torch.randn(input_dim, latent_dim, latent_dim, dtype=torch.complex64))
        )
        self.S = nn.Parameter(
            0.1 * (torch.randn(output_dim, latent_dim, latent_dim, dtype=torch.complex64))
        )

        self.bias = nn.Parameter(torch.zeros(output_dim, dtype=torch.float32))

    @staticmethod
    def hermitian(A):
        """
        Make A Hermitian:
            (A + A^†) / 2
        Supports batch: (..., n, n)
        """
        return 0.5 * (A + A.conj().transpose(-1, -2))

    def forward(self, x):
        """
        x: real input tensor of shape (batch, input_dim).
        Returns:
            y: real tensor of shape (batch, output_dim).
        """
        batch_size = x.shape[0]
        x_real = x.to(torch.float32)

        # Hermitian primary params
        H0_herm = self.hermitian(self.H0)              # (n, n)
        H_inputs_herm = self.hermitian(self.H_inputs)  # (p, n, n)

        # Build H(x)
        H = H0_herm.unsqueeze(0).expand(batch_size, -1, -1).clone()  # (b, n, n)
        H = H.to(torch.complex64)
        for l in range(self.input_dim):
            coef = x_real[:, l].view(-1, 1, 1).to(torch.complex64)
            H = H + coef * H_inputs_herm[l]

        H = self.hermitian(H)

        # Eigen-decomposition
        eigvals, eigvecs = torch.linalg.eigh(H)  # eigvecs: (b, n, n)

        # Take first r eigenvectors
        v = eigvecs[:, :, :self.r]  # (b, n, r)

        # Hermitian secondary matrices
        S_herm = self.hermitian(self.S)  # (out, n, n)

        # Output
        y = torch.zeros(batch_size, self.output_dim, dtype=torch.float32, device=x.device)
        for k in range(self.output_dim):
            S_k = S_herm[k]  # (n, n)
            val_k = torch.zeros(batch_size, dtype=torch.complex64, device=x.device)
            for i in range(self.r):
                v_i = v[:, :, i]  # (b, n)
                tmp = torch.einsum("bi,ij,bj->b", v_i.conj(), S_k, v_i)  # (b,)
                val_k += tmp
            y[:, k] = val_k.real + self.bias[k]

        return y


# ============================================================
# 4. Training + diagnostics + trajectory emulation
# ============================================================

def main():
    # ----------------------------
    # Generate sequence dataset
    # ----------------------------
    print("Generating sequence dataset from exact Lipkin dynamics...")
    N = 10         # number of particles (dimension = N+1)
    t_max = 10.0
    n_steps = 200
    epsilon = 1.0

    X, Y = generate_sequence_dataset(
        num_trajectories=60,
        N=N,
        n_steps=n_steps,
        t_max=t_max,
        epsilon=epsilon,
        chi_min=0.1,
        chi_max=1.5,
        subsample=1,
        seed=1234
    )

    N_pairs, input_dim = X.shape
    output_dim = Y.shape[1]
    print(f"Number of one-step pairs: {N_pairs}")
    print(f"Input dimension: {input_dim}, Output dimension: {output_dim}")

    # Train/validation split
    frac_train = 0.8
    N_train = int(frac_train * N_pairs)
    perm = np.random.permutation(N_pairs)
    train_idx = perm[:N_train]
    val_idx = perm[N_train:]

    X_train = torch.tensor(X[train_idx], dtype=torch.float32)
    Y_train = torch.tensor(Y[train_idx], dtype=torch.float32)
    X_val = torch.tensor(X[val_idx], dtype=torch.float32)
    Y_val = torch.tensor(Y[val_idx], dtype=torch.float32)

    # ----------------------------
    # Instantiate PMM sequence model
    # ----------------------------
    model = UnitaryPMMSequence(
        input_dim=input_dim,
        output_dim=output_dim,
        latent_dim=8,
        r=3
    )

    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.MSELoss()

    # ----------------------------
    # Training loop
    # ----------------------------
    print("Training PMM sequence model...")
    num_epochs = 300
    batch_size = 128

    N_train_torch = X_train.shape[0]
    for epoch in range(1, num_epochs + 1):
        model.train()
        perm_torch = torch.randperm(N_train_torch)
        total_loss = 0.0
        n_batches = 0

        for i in range(0, N_train_torch, batch_size):
            idx = perm_torch[i:i+batch_size]
            x_batch = X_train[idx]
            y_batch = Y_train[idx]

            optimizer.zero_grad()
            y_pred = model(x_batch)
            loss = loss_fn(y_pred, y_batch)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            n_batches += 1

        train_loss = total_loss / max(1, n_batches)

        # Validation
        model.eval()
        with torch.no_grad():
            y_val_pred = model(X_val)
            val_loss = loss_fn(y_val_pred, Y_val).item()

        if epoch % 50 == 0 or epoch == 1:
            print(f"Epoch {epoch:4d}: "
                  f"train MSE = {train_loss:.4e}, "
                  f"val MSE = {val_loss:.4e}")

    print("Training finished.")

    # ======================================================
    # 5. Energy conservation diagnostics for a sample traj
    # ======================================================
    print("\nEnergy conservation diagnostics...")

    rng = np.random.default_rng(2025)
    chi_test = rng.uniform(0.1, 1.5)
    theta0 = np.arccos(2 * rng.random() - 1.0)
    phi0 = 2 * np.pi * rng.random()

    model_q = LipkinQuantum(N=N, epsilon=epsilon, chi=chi_test)
    psi0 = model_q.coherent_state(theta0, phi0)
    times, psis = model_q.evolve_rk4(psi0, t_max=t_max, n_steps=n_steps)

    E = np.array([model_q.expectation_energy(psis[i]) for i in range(n_steps + 1)])
    E0 = E[0]
    rel_dev = (E - E0) / (abs(E0) + 1e-14)

    fig_E, axes_E = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    axes_E[0].plot(times, E)
    axes_E[0].set_ylabel("Energy E(t)")
    axes_E[0].grid(True)

    axes_E[1].plot(times, rel_dev)
    axes_E[1].set_xlabel("time t")
    axes_E[1].set_ylabel("rel. dev. (E(t)-E(0))/|E(0)|")
    axes_E[1].grid(True)

    fig_E.suptitle("Energy conservation in finite-N Lipkin dynamics")
    plt.tight_layout()

    # ======================================================
    # 6. Trajectory emulation by iterated PMM step
    # ======================================================
    print("\nEmulating a trajectory via iterated PMM step...")

    # True Bloch vectors from exact dynamics
    S_true = np.zeros((n_steps + 1, 3), dtype=np.float64)
    for i in range(n_steps + 1):
        S_true[i] = model_q.expectation_J(psis[i])

    dt = t_max / n_steps
    # Initialize PMM trajectory with the exact initial Bloch vector
    S_pmm = np.zeros_like(S_true)
    S_pmm[0] = S_true[0]

    # Iterate PMM step
    model.eval()
    with torch.no_grad():
        for n in range(n_steps):
            s_n = S_pmm[n]
            x_n = np.array([
                s_n[0],
                s_n[1],
                s_n[2],
                epsilon,
                chi_test,
                dt / t_max
            ], dtype=np.float32)
            x_t = torch.tensor(x_n, dtype=torch.float32).unsqueeze(0)  # (1, input_dim)
            y_pred = model(x_t).cpu().numpy()[0]  # (3,)
            S_pmm[n+1] = y_pred

    # Plot Bloch components: exact vs PMM emulation
    fig_S, axes_S = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    labels = [r"$s_x$", r"$s_y$", r"$s_z$"]
    for i, ax in enumerate(axes_S):
        ax.plot(times, S_true[:, i], label="Exact quantum", linewidth=2)
        ax.plot(times, S_pmm[:, i], "--", label="PMM sequence", linewidth=1.5)
        ax.set_ylabel(labels[i])
        ax.grid(True)
        if i == 0:
            ax.legend()
    axes_S[-1].set_xlabel("time t")
    fig_S.suptitle("Finite-N Lipkin dynamics: exact vs PMM sequence emulation")
    plt.tight_layout()

    plt.show()

    # Print a few numerical comparisons
    print("\nSample numerical comparisons (first 5 time steps):")
    for k in range(5):
        print(f"t = {times[k]:6.3f}: "
              f"s_true = {S_true[k]}, "
              f"s_pmm = {S_pmm[k]}")


if __name__ == "__main__":
    main()
