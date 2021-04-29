import os
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from sksparse.cholmod import cholesky
import sys

build_directory = "../"

mesh_name = sys.argv[1]  # 'tube_1x1_h0.08'
out_dir = build_directory + 'basis/'
matrix_dir = build_directory + 'matrices/'
tol_velocity = 1e-4
tol_pressure = 1e-5

generate_velocity = True
generate_pressure = True
generate_primal_supremizers = True
generate_dual_supremizers = True
generate_supremizers = generate_primal_supremizers or generate_dual_supremizers

single_block_flag = True  # True if the simulation is performed on a single block only
weak_dirichlet_inflow = True  # True if the BC at the inlet is a Dirichlet BC with weak imposition via multipliers

# this is the maximum number of snapshots we want to process
n_snapshots = 200

# if True, we normalize with respect to custom norms; else we use standard l2-norm
use_custom_norms = True


def create_dir(name):
    try:
        os.mkdir(name)
    except FileExistsError:
        pass


def generate_basis(index, nsnaps, norm_matrix, tol):
    print("===== COMPUTING BASIS FIELD " + str(index) + " ====", flush=True)

    if use_custom_norms:
        factor = cholesky(norm_matrix)

    snap = False
    count = 0
    print("reading snapshots ...", flush=True)
    for i in range(nsnaps):
        fname = build_directory + 'snapshots/param' + str(i) + '/' + mesh_name + '/field' + str(index) + '.snap'
        if os.path.isfile(fname):
            count = count + 1
            print('\t snapshotfile: ' + fname, flush=True)
            cur_data = np.genfromtxt(fname, delimiter=',')
            if snap is False:
                snap = cur_data.T
            else:
                snap = np.concatenate((snap, cur_data.T), axis=1)

    print("Number of read snapshots = " + str(count), flush=False)
    print("\n")
    print('Size of snapshots matrix = ' + str(snap.shape), flush=True)

    if use_custom_norms:
        snap = factor.L().transpose().dot(factor.apply_P(snap))
    U, S, V = np.linalg.svd(snap, full_matrices=False)
    if use_custom_norms:
        U = factor.apply_Pt(factor.solve_Lt(U, use_LDLt_decomposition=False))

    total_energy = np.sum(np.square(S))
    print('Total energy = ' + str(total_energy), flush=True)

    print("\n")
    partial_sum = 0
    Nu = 0
    while tol * tol < 1.0 - partial_sum / total_energy:
        partial_sum = partial_sum + S[Nu] * S[Nu]
        msg = 'index = ' + str(Nu)
        msg += ' partial sum = ' + str(partial_sum)
        msg += ' remainder = ' + str(1.0 - partial_sum / total_energy)
        print(msg)
        Nu = Nu + 1

    U = U[:, :Nu]

    print("dumping to file ...", flush=True)

    create_dir(out_dir)
    create_dir(out_dir + '/' + mesh_name)

    U[np.abs(U) < 1e-15] = 0
    np.savetxt(out_dir + '/' + mesh_name + '/field' + str(index) + '.basis', U.T,
               fmt='%.16g', delimiter=',')
    np.savetxt(out_dir + '/' + mesh_name + '/svd' + str(index) + '.txt', S[:, np.newaxis],
               fmt='%.16g', delimiter=',')
    print("done\n", flush=True)

    return U


def read_matrix(name, matrixtype=csr_matrix, shapem=None):
    print(f"Reading matrix {name}\n")
    ii, jj, value = np.loadtxt(name).T
    # indices are relative to matlab so we subtract 1
    return matrixtype((value, (ii.astype(int) - 1, jj.astype(int) - 1)), shape=shapem)


def mydot(vec1, vec2, norm_matrix):
    if use_custom_norms:
        aux = norm_matrix.dot(vec2)
        return vec1.dot(aux)
    return vec1.dot(vec2)


def mynorm(vec, norm_matrix):
    return np.sqrt(mydot(vec, vec, norm_matrix))


norm_velocity = read_matrix(matrix_dir + '/' + mesh_name + '/norm0.m') if generate_supremizers else None
norm_pressure = read_matrix(matrix_dir + '/' + mesh_name + '/norm1.m') if generate_supremizers else None

# find Dirichlet nodes
if generate_supremizers:
    D = norm_velocity.diagonal()
    dir_indices = np.abs(D - 1) < 1e-16

norm_velocity_nobcs = read_matrix(matrix_dir + '/' + mesh_name + '/norm0_nobcs.m') if use_custom_norms else None
norm_pressure_nobcs = read_matrix(matrix_dir + '/' + mesh_name + '/norm1_nobcs.m') if use_custom_norms else None

if generate_velocity:
    U = generate_basis(0, n_snapshots, norm_velocity_nobcs, tol_velocity)
else:
    print("===== IMPORTING BASIS FIELD 0 ====", flush=True)
    U = np.genfromtxt(out_dir + '/' + mesh_name + '/field0.basis', delimiter=',').T
if generate_pressure:
    P = generate_basis(1, n_snapshots, norm_pressure_nobcs, tol_pressure)
else:
    print("===== IMPORTING BASIS FIELD 1 ====", flush=True)
    P = np.genfromtxt(out_dir + '/' + mesh_name + '/field1.basis', delimiter=',').T

if generate_primal_supremizers:
    print("===== COMPUTING PRIMAL SUPREMIZERS ====", flush=True)

    constraint_matrix = read_matrix(matrix_dir + '/' + mesh_name + '/primalConstraint.m', csr_matrix)
    constraint_matrix[dir_indices] = 0
    rhs = constraint_matrix * P
    supr_primal = linalg.spsolve(norm_velocity, rhs)

    minnorm = 1e16
    for i in range(supr_primal.shape[1]):
        print('Normalizing primal supremizer ' + str(i), flush=True)
        # # normalize
        for col in U.T:
            supr_primal[:, i] = supr_primal[:, i] - \
                                mydot(supr_primal[:, i], col, norm_velocity_nobcs) / \
                                mynorm(col, norm_velocity_nobcs) * col

        for j in range(i):
            supr_primal[:, i] = supr_primal[:, i] - \
                                mydot(supr_primal[:, i], supr_primal[:, j], norm_velocity_nobcs) / \
                                mynorm(supr_primal[:, j], norm_velocity_nobcs) * supr_primal[:, j]

        nnorm = mynorm(supr_primal[:, i], norm_velocity_nobcs)
        minnorm = min(nnorm, minnorm)
        print('\tnorm supremizers = ' + str(nnorm))
        supr_primal[:, i] = supr_primal[:, i] / nnorm
    print('Min norm = ' + str(minnorm))

    supr_primal[np.abs(supr_primal) < 1e-15] = 0
    np.savetxt(out_dir + '/' + mesh_name + '/primal_supremizers_0_1.basis', supr_primal.T, fmt='%.16g', delimiter=',')
elif generate_dual_supremizers:
    try:
        print("===== IMPORTING PRIMAL SUPREMIZERS ====", flush=True)
        supr_primal = np.genfromtxt(out_dir + '/' + mesh_name + '/primal_supremizers_0_1.basis', delimiter=',').T
    except (FileNotFoundError, OSError):
        print("===== PRIMAL SUPREMIZERS IMPORT FAILED ====", flush=True)
        supr_primal = np.zeros(0)

if generate_dual_supremizers:
    print("===== COMPUTING DUAL SUPREMIZERS ====", flush=True)

    constraint_matrices = []

    if weak_dirichlet_inflow:
        constraint_matrices.append(read_matrix(matrix_dir + '/' + mesh_name + '/dualConstraint1.m', csc_matrix))  # inflow

    if not single_block_flag:
        constraint_matrices.append(read_matrix(matrix_dir + '/' + mesh_name + '/dualConstraint2.m', csc_matrix))  # outflow 1
        if mesh_name == 'bif_sym_alpha50_0.10':
            constraint_matrices.append(read_matrix(matrix_dir + '/' + mesh_name + '/dualConstraint3.m', csc_matrix))  # outflow 2

    for constraint_matrix in constraint_matrices:
        pad = csc_matrix((U.shape[0] - constraint_matrix.shape[0], constraint_matrix.shape[1]))
        constraint_matrix = scipy.sparse.vstack([constraint_matrix, pad])

    global_constraint = csr_matrix(scipy.sparse.hstack([constraint_matrix for constraint_matrix in constraint_matrices]))
    global_constraint[dir_indices] = 0

    print(global_constraint.shape)
    supr_dual = linalg.spsolve(norm_velocity, global_constraint).toarray()
    minnorm = 1e16
    for i in range(supr_dual.shape[1]):
        print('Normalizing dual supremizer ' + str(i), flush=True)

        for col in U.T:
            supr_dual[:, i] = supr_dual[:, i] - \
                              mydot(supr_dual[:, i], col, norm_velocity_nobcs) / \
                              mynorm(col, norm_velocity_nobcs) * col

        for col in supr_primal.T:
            supr_dual[:, i] = supr_dual[:, i] - \
                              mydot(supr_dual[:, i], col, norm_velocity_nobcs) / \
                              mynorm(col, norm_velocity_nobcs) * col

        for j in range(i):
            supr_dual[:, i] = supr_dual[:, i] - \
                              mydot(supr_dual[:, i], supr_dual[:, j], norm_velocity_nobcs) / \
                              mynorm(supr_dual[:, j], norm_velocity_nobcs) * supr_dual[:, j]

        nnorm = mynorm(supr_dual[:, i], norm_velocity_nobcs)
        minnorm = min(nnorm, minnorm)
        print('\tnorm supremizers = ', str(nnorm))
        supr_dual[:, i] = supr_dual[:, i] / nnorm
    print('Min norm = ' + str(minnorm))

    supr_dual[np.abs(supr_dual) < 1e-15] = 0
    np.savetxt(out_dir + '/' + mesh_name + '/dual_supremizers0.basis', supr_dual.T, fmt='%.16g', delimiter=',')
