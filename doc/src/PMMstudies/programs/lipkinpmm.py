"""
Uses a more realistic Lipkin–Meshkov–Glick mean-field (TDHF) Hamiltonian h[\mathbf{s}] = \varepsilon J_z - \frac{\chi}{2} J_x^2 \Rightarrow classical spin equations \dot{\mathbf{s}} = \mathbf{s} \times \mathbf{B}_\mathrm{eff} with \mathbf{B}_\mathrm{eff} = (-\chi s_x,\, 0,\, \varepsilon).
Implements TDHF-like evolution as Bloch vector dynamics.
Implements a more faithful PMM: 
Primary matrix is a Hermitian matrix depending affinely on the inputs.
It is interpreted as an effective generator of a unitary via U = e^{-iH}, though we only need its eigenvectors (same as the unitary’s).
Outputs are expectation values of Hermitian secondary matrices in a few low-lying eigenvectors.

Trains the PMM to emulate the TDHF dynamics.
Plots TDHF vs PMM predictions for one test trajectory (for s_x, s_y, s_z vs time).

TDHF-like evolution + Parametric Matrix Model (PMM) emulator
for the Lipkin–Meshkov–Glick (Lipkin) model.

Features:
  * More realistic Lipkin mean-field Hamiltonian:
        H/N = epsilon * s_z - (chi/2) * s_x^2
    leading to Bloch equations: ds/dt = s x B_eff,
        B_eff = ( -chi * s_x, 0, epsilon )

  * TDHF evolution implemented as classical spin (Bloch vector) dynamics.
  * PMM with complex Hermitian primary and secondary matrices, interpreted
    as a unitary-based PMM (H is generator of a unitary U = exp(-i H)).
  * Training of the PMM to emulate TDHF trajectories.
  * Plotting of TDHF vs PMM predictions for one test trajectory.


"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt


# ============================================================
# 1. Lipkin TDHF-like evolution (classical spin / Bloch vector)
# ============================================================

class LipkinTDHFRealistic:
    """
    TDHF-like mean-field dynamics for the Lipkin–Meshkov–Glick model
    in the large-N (classical spin) limit.

    We use a normalized Bloch vector s = (sx, sy, sz) with |s| ~ 1
    and a Hamiltonian density:
        h(s) = epsilon * s_z - (chi/2) * s_x^2

    The effective field is:
        B_eff = ( -chi * s_x, 0, epsilon )

    Equations of motion (classical spin / TDHF mean-field):
        ds/dt = s x B_eff
    """

    def __init__(self, epsilon=1.0, chi=1.0):
        self.epsilon = epsilon
        self.chi = chi

    def rhs(self, s):
        """
        Right-hand side ds/dt = s x B_eff
        s: array-like shape (3,)
        Returns:
            dsdt: numpy array shape (3,)
        """
        sx, sy, sz = s
        Bx = -self.chi * sx
        By = 0.0
        Bz = self.epsilon

        # cross product s x B
        dsdt = np.array([
            sy * Bz - sz * By,
            sz * Bx - sx * Bz,
            sx * By - sy * Bx
        ], dtype=np.float64)

        return dsdt

    def evolve_rk4(self, s0, t_max, n_steps):
        """
        Evolve Bloch vector using 4th-order Runge-Kutta.

        Parameters:
            s0     : initial Bloch vector (3,)
            t_max  : final time
            n_steps: number of time steps

        Returns:
            times: array of times, shape (n_steps+1,)
            S    : array of Bloch vectors, shape (n_steps+1, 3)
        """
        dt = t_max / n_steps
        times = np.linspace(0.0, t_max, n_steps + 1)
        S = np.zeros((n_steps + 1, 3), dtype=np.float64)
        S[0] = s0

        s = s0.copy()
        for n in range(n_steps):
            k1 = self.rhs(s)
            k2 = self.rhs(s + 0.5 * dt * k1)
            k3 = self.rhs(s + 0.5 * dt * k2)
            k4 = self.rhs(s + dt * k3)

            s = s + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

            # Optionally renormalize to |s| ~ 1 to avoid drift
            norm = np.linalg.norm(s)
            if norm > 1e-12:
                s = s / norm

            S[n+1] = s

        return times, S


# ============================================================
# 2. Dataset generation from TDHF trajectories
# ============================================================

def generate_tdhf_dataset(num_trajectories=40,
                          n_steps=200,
                          t_max=10.0,
                          epsilon=1.0,
                          chi_min=0.1,
                          chi_max=1.5,
                          subsample=4,
                          seed=1234):
    """
    Generate a dataset from TDHF-like trajectories of the Lipkin model.

    For each trajectory:
      - Sample chi uniformly in [chi_min, chi_max].
      - Sample initial state on the Bloch sphere (theta, phi).
      - Evolve LipkinTDHFRealistic.
      - For selected times, store:
            inputs  x = [t/t_max, epsilon, chi,
                        cos(theta0), sin(theta0),
                        cos(phi0),   sin(phi0)]
            outputs y = [sx(t), sy(t), sz(t)]

    Returns:
        X: numpy array (N_samples, input_dim)
        y: numpy array (N_samples, 3)
    """
    rng = np.random.default_rng(seed)

    X_list = []
    y_list = []

    for traj in range(num_trajectories):
        # Sample chi and initial angles
        chi = rng.uniform(chi_min, chi_max)
        theta0 = np.arccos(2 * rng.random() - 1.0)
        phi0 = 2 * np.pi * rng.random()

        s0 = np.array([
            np.sin(theta0) * np.cos(phi0),
            np.sin(theta0) * np.sin(phi0),
            np.cos(theta0)
        ], dtype=np.float64)

        model = LipkinTDHFRealistic(epsilon=epsilon, chi=chi)
        times, S = model.evolve_rk4(s0, t_max=t_max, n_steps=n_steps)

        for n in range(0, n_steps + 1, subsample):
            t = times[n]
            s = S[n]

            x = np.array([
                t / t_max,
                epsilon,
                chi,
                np.cos(theta0),
                np.sin(theta0),
                np.cos(phi0),
                np.sin(phi0)
            ], dtype=np.float64)

            X_list.append(x)
            y_list.append(s)

    X = np.stack(X_list, axis=0)
    y = np.stack(y_list, axis=0)
    return X, y


# ============================================================
# 3. Complex Hermitian PMM (unitary-based)
# ============================================================

class UnitaryPMM(nn.Module):
    """
    More faithful PMM:
      - Primary matrices are complex Hermitian:
            H(x) = H0 + sum_l x_l H_l
      - Interpreted as generator of a unitary U(x) = exp(-i H(x)).
        (We only need eigenvectors of H, which are also those of U.)
      - Secondary matrices S_k are complex Hermitian.
      - Outputs are real expectation values:
            y_k(x) = bias_k + sum_{i=1}^r Re( v_i^\dagger S_k v_i )

    where v_i are the r eigenvectors corresponding to the lowest eigenvalues
    (or first r eigenvalues) of H(x).
    """

    def __init__(self, input_dim, output_dim, latent_dim=6, r=3):
        super().__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.latent_dim = latent_dim
        self.r = r

        # Primary matrix parameters: complex
        # H0 and H_inputs[l] are unconstrained complex matrices;
        # we will symmetrize them to get Hermitian matrices.
        self.H0 = nn.Parameter(
            0.1 * (torch.randn(latent_dim, latent_dim) +
                   1j * torch.randn(latent_dim, latent_dim))
        )
        self.H_inputs = nn.Parameter(
            0.1 * (torch.randn(input_dim, latent_dim, latent_dim) +
                   1j * torch.randn(input_dim, latent_dim, latent_dim))
        )

        # Secondary matrices: complex, will be symmetrized
        self.S = nn.Parameter(
            0.1 * (torch.randn(output_dim, latent_dim, latent_dim) +
                   1j * torch.randn(output_dim, latent_dim, latent_dim))
        )

        # Real bias term
        self.bias = nn.Parameter(torch.zeros(output_dim))

    @staticmethod
    def hermitian(A):
        """
        Make A (possibly batch of matrices) Hermitian:
            (A + A^†) / 2
        """
        return 0.5 * (A + A.conj().transpose(-1, -2))

    def forward(self, x):
        """
        x: (batch, input_dim), real tensor.
        Returns:
            y: (batch, output_dim), real tensor.
        """
        batch_size = x.shape[0]

        # Build Hermitian H0 and H_inputs
        H0_herm = self.hermitian(self.H0)  # (n, n)
        H_inputs_herm = self.hermitian(self.H_inputs)  # (p, n, n)

        # Build H(x) for each batch element
        H = H0_herm.unsqueeze(0).expand(batch_size, -1, -1).clone()  # (b, n, n)
        x_c = x.to(dtype=torch.float32)  # real; we'll cast scalars to complex
        for l in range(self.input_dim):
            coef = x_c[:, l].view(-1, 1, 1).to(H.dtype)
            H = H + coef * H_inputs_herm[l]

        H = self.hermitian(H)

        # Eigen-decompose H(x): H v = lambda v
        eigvals, eigvecs = torch.linalg.eigh(H)  # eigvecs: (b, n, n)

        # Take the first r eigenvectors (smallest eigenvalues)
        v = eigvecs[:, :, :self.r]  # (b, n, r)

        # Hermitian secondary matrices
        S_herm = self.hermitian(self.S)  # (out, n, n)

        # Compute outputs
        y = torch.zeros(batch_size, self.output_dim, dtype=torch.float32, device=x.device)

        # For each output component k
        for k in range(self.output_dim):
            S_k = S_herm[k]  # (n, n)
            val_k = torch.zeros(batch_size, dtype=torch.complex64, device=x.device)
            # Sum over r eigenvectors
            for i in range(self.r):
                v_i = v[:, :, i]  # (b, n)
                # v_i^\dagger S_k v_i -> shape (b,)
                tmp = torch.einsum("bi,ij,bj->b", v_i.conj(), S_k, v_i)
                val_k += tmp
            # Take real part and add bias
            y[:, k] = val_k.real + self.bias[k]

        return y


# ============================================================
# 4. Training + plotting
# ============================================================

def main():
    # ----------------------------
    # Generate TDHF dataset
    # ----------------------------
    print("Generating TDHF dataset...")
    X, y = generate_tdhf_dataset(
        num_trajectories=60,
        n_steps=200,
        t_max=10.0,
        epsilon=1.0,
        chi_min=0.1,
        chi_max=1.5,
        subsample=3,
        seed=1234
    )

    N_samples, input_dim = X.shape
    output_dim = y.shape[1]
    print(f"Dataset size: {N_samples}")
    print(f"Input dimension: {input_dim}, Output dimension: {output_dim}")

    # Train/validation split
    frac_train = 0.8
    N_train = int(frac_train * N_samples)
    perm_all = np.random.permutation(N_samples)
    train_idx = perm_all[:N_train]
    val_idx = perm_all[N_train:]

    X_train = torch.tensor(X[train_idx], dtype=torch.float32)
    y_train = torch.tensor(y[train_idx], dtype=torch.float32)
    X_val = torch.tensor(X[val_idx], dtype=torch.float32)
    y_val = torch.tensor(y[val_idx], dtype=torch.float32)

    # ----------------------------
    # Instantiate PMM model
    # ----------------------------
    latent_dim = 8
    r = 3
    model = UnitaryPMM(
        input_dim=input_dim,
        output_dim=output_dim,
        latent_dim=latent_dim,
        r=r
    )

    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.MSELoss()

    # ----------------------------
    # Training loop
    # ----------------------------
    num_epochs = 400
    batch_size = 128

    print("Starting training...")
    N_train_torch = X_train.shape[0]
    for epoch in range(1, num_epochs + 1):
        model.train()
        perm = torch.randperm(N_train_torch)
        total_loss = 0.0
        n_batches = 0

        for i in range(0, N_train_torch, batch_size):
            idx = perm[i:i+batch_size]
            x_batch = X_train[idx]
            y_batch = y_train[idx]

            optimizer.zero_grad()
            y_pred = model(x_batch)
            loss = loss_fn(y_pred, y_batch)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            n_batches += 1

        train_loss = total_loss / max(1, n_batches)

        model.eval()
        with torch.no_grad():
            y_val_pred = model(X_val)
            val_loss = loss_fn(y_val_pred, y_val).item()

        if epoch % 50 == 0 or epoch == 1:
            print(f"Epoch {epoch:4d}: "
                  f"train MSE = {train_loss:.4e}, "
                  f"val MSE = {val_loss:.4e}")

    print("Training finished.")

    # =======================================================
    # 5. Compare TDHF vs PMM on an independent trajectory
    # =======================================================
    print("\nGenerating independent test trajectory for plotting...")

    # Independent chi and initial state
    rng = np.random.default_rng(2025)
    epsilon = 1.0
    chi = rng.uniform(0.1, 1.5)
    theta0 = np.arccos(2 * rng.random() - 1.0)
    phi0 = 2 * np.pi * rng.random()
    s0 = np.array([
        np.sin(theta0) * np.cos(phi0),
        np.sin(theta0) * np.sin(phi0),
        np.cos(theta0)
    ], dtype=np.float64)

    t_max = 10.0
    n_steps = 200
    lipkin_test = LipkinTDHFRealistic(epsilon=epsilon, chi=chi)
    times, S_true = lipkin_test.evolve_rk4(s0, t_max=t_max, n_steps=n_steps)

    # Build inputs for PMM for this trajectory
    X_test = []
    for t in times:
        x = np.array([
            t / t_max,
            epsilon,
            chi,
            np.cos(theta0),
            np.sin(theta0),
            np.cos(phi0),
            np.sin(phi0)
        ], dtype=np.float64)
        X_test.append(x)
    X_test = np.stack(X_test, axis=0)
    X_test_torch = torch.tensor(X_test, dtype=torch.float32)

    model.eval()
    with torch.no_grad():
        S_pred = model(X_test_torch).cpu().numpy()  # shape (n_steps+1, 3)

    # ----------------------------
    # Plot results
    # ----------------------------
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

    labels = [r"$s_x$", r"$s_y$", r"$s_z$"]
    for i, ax in enumerate(axes):
        ax.plot(times, S_true[:, i], label="TDHF (Lipkin)")
        ax.plot(times, S_pred[:, i], linestyle="--", label="PMM prediction")
        ax.set_ylabel(labels[i])
        ax.grid(True)
        if i == 0:
            ax.legend()

    axes[-1].set_xlabel("time t")

    fig.suptitle("Lipkin TDHF vs PMM emulation (Bloch components)")
    plt.tight_layout()
    plt.show()

    # Also print a few numerical comparisons
    print("\nExample numerical comparison (first 5 times):")
    for k in range(5):
        print(f"t = {times[k]:6.3f}: "
              f"s_true = {S_true[k]}, "
              f"s_pred = {S_pred[k]}")


if __name__ == "__main__":
    main()
