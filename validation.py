# validation.py
import numpy as np


from code_to_be_implemented.elements import tetra4_volume, tetra4_strain_stress, tetra4_K, hooke_iso_3d
from code_to_be_implemented.assembly import assembly
from code_to_be_implemented.boundcond import boundcond
from inputs import input_single_element
from inputs import input_cube

node_single = input_single_element.node
mater_single = input_single_element.mater[0]

# SINGLE ELEMENT
# 1.1. Volume
vol = tetra4_volume(node_single)
expected_vol = 1.0/6.0
if np.isclose(vol, expected_vol):
    res = "OK"
else:
    res = "FAIL"
print(f"SINGLE ELEMENT - Volume ; EXPECTED : {expected_vol:.3f} ; CALCULATED : {vol:.3f} -> {res}")

# 1.2. Patch tests
# 1.2.1. Zero
u_zero = np.zeros(12)
strain_zero, stress_zero = tetra4_strain_stress(node_single, mater_single, u_zero)
max_strain = np.max(np.abs(strain_zero))
if np.allclose(strain_zero, 0) and np.allclose(stress_zero, 0):
    res = "OK"
else:
    res = "FAIL"
print(f"SINGLE ELEMENT - Patch test (zero) ; EXPECTED : 0.0 ; CALCULATED : maxx={max_strain:.3f} -> {res}")

# 1.2.2 Constant strain tests
strain_labels = ["exx", "eyy", "ezz", "gxy", "gyz", "gxz"]

for i in range(6):
    label = strain_labels[i]

    u = np.zeros(12)
    for n in range(4):
        x = node_single[n,0]
        y = node_single[n,1]
        z = node_single[n,2]

        if i == 0:              # exx
            u[3*n] = x
        elif i == 1:            # eyy
            u[3*n+1] = y
        elif i == 2:            # ezz
            u[3*n+2] = z
        elif i == 3:            # gxy
            u[3*n] = 0.5*y
            u[3*n+1] = 0.5*x
        elif i == 4:            # gyz
            u[3*n+1] = 0.5*z
            u[3*n+2] = 0.5*y
        elif i == 5:            # gxz
            u[3*n] = 0.5*z
            u[3*n+2] = 0.5*x

    strain, stress = tetra4_strain_stress(node_single, mater_single, u)

    expected_strain = np.zeros(6)
    expected_strain[i] = 1.0

    H_mat = hooke_iso_3d(mater_single[0], mater_single[1])
    expected_stress = H_mat @ expected_strain

    err_strain = np.max(np.abs(strain - expected_strain))
    err_stress = np.max(np.abs(stress - expected_stress))

    if np.allclose(strain, expected_strain) and np.allclose(stress, expected_stress):
        res = "OK"
    else:
        res = "FAIL"
    print(f"SINGLE ELEMENT - Patch test ({label}) ; max_error_strain={err_strain:.3e}, max_error_stress={err_stress:.3e} -> {res}")

# 1.3. Rigid body (Translations)
axes = ["X", "Y", "Z"]
for i in range(3):
    axis = axes[i]
    u = np.zeros(12)
    for n in range(4):
        u[3*n + i] = 1.0

    strain, stress = tetra4_strain_stress(node_single, mater_single, u)
    max_val = np.max(np.abs(strain))

    if np.allclose(strain, 0) and np.allclose(stress, 0):
        res = "OK"
    else:
        res = "FAIL"
    print(f"SINGLE ELEMENT - Rigid body (Trans {axis}) ; EXPECTED : 0.0 ; CALCULATED : max={max_val:.3f} -> {res}")

# 1.4. Rigid body (Rotations)
for i in range(3):
    axis = axes[i]
    u = np.zeros(12)
    for n in range(4):
        x = node_single[n,0]
        y = node_single[n,1]
        z = node_single[n,2]

        if i == 0:              # Rot X
            u[3*n+1] = -z
            u[3*n+2] = y
        elif i == 1:            # Rot Y
            u[3*n] = z
            u[3*n+2] = -x
        elif i == 2:            # Rot Z
            u[3*n] = -y
            u[3*n+1] = x

    strain, stress = tetra4_strain_stress(node_single, mater_single, u)
    max_val = np.max(np.abs(strain))

    if np.allclose(strain, 0) and np.allclose(stress, 0):
        res = "OK"
    else:
        res = "FAIL"
    print(f"SINGLE ELEMENT - Rigid body (Rot {axis}) ; EXPECTED : 0.0 ; CALCULATED : max={max_val:.3f} -> {res}")

# 1.5. Poisson ratio
mater_nu0 = np.array([1000.0, 0.0])
u_exx = np.zeros(12)
for n in range(4):
    u_exx[3*n] = node_single[n,0]

strain_nu0, stress_nu0 = tetra4_strain_stress(node_single, mater_nu0, u_exx)
if np.isclose(stress_nu0[1], 0.0) and np.isclose(stress_nu0[2], 0.0):
    res = "OK"
else:
    res = "FAIL"
print(f"SINGLE ELEMENT - Poisson ratio (nu=0) ; EXPECTED : 0.0 ; CALCULATED : {stress_nu0[1]:.3f} -> {res}")

Ksys_cube = assembly(input_cube.node, input_cube.elem, input_cube.eltp, input_cube.mater)
ndof_cube = Ksys_cube.shape[0]

# MULTI ELEMENT
# 2.1. Patch test
for i in range(6):
    label = strain_labels[i]

    u_global = np.zeros(ndof_cube)
    for n in range(input_cube.node.shape[0]):
        x = input_cube.node[n,0]
        y = input_cube.node[n,1]
        z = input_cube.node[n,2]

        if i == 0:                  # exx
            u_global[3*n] = x
        elif i == 1:                # eyy
            u_global[3*n+1] = y
        elif i == 2:                # ezz
            u_global[3*n+2] = z
        elif i == 3:                # gxy
            u_global[3*n] = 0.5*y
            u_global[3*n+1] = 0.5*x
        elif i == 4:                # gyz
            u_global[3*n+1] = 0.5*z
            u_global[3*n+2] = 0.5*y
        elif i == 5:                # gxz
            u_global[3*n] = 0.5*z
            u_global[3*n+2] = 0.5*x

    max_err_strain = 0.0
    max_err_stress = 0.0
    all_ok = True

    for e in range(input_cube.elem.shape[0]):
        node_ids = input_cube.elem[e, 2:6] - 1
        node_e = input_cube.node[node_ids]
        u_e = np.zeros(12)
        for j in range(4):
            u_e[3*j:3*j+3] = u_global[3*node_ids[j]:3*node_ids[j]+3]

        strain, stress = tetra4_strain_stress(node_e, input_cube.mater[0], u_e)

        expected_strain = np.zeros(6)
        expected_strain[i] = 1.0

        H_mat = hooke_iso_3d(input_cube.mater[0,0], input_cube.mater[0,1])
        expected_stress = H_mat @ expected_strain

        err_strain = np.max(np.abs(strain - expected_strain))
        if err_strain > max_err_strain:
            max_err_strain = err_strain

        err_stress = np.max(np.abs(stress - expected_stress))
        if err_stress > max_err_stress:
            max_err_stress = err_stress

        if not np.allclose(strain, expected_strain) and np.allclose(stress, expected_stress):
            all_ok = False

    if all_ok:
        res = "OK"
    else:
        res = "FAIL"
    print(f"MULTI ELEMENT - Patch test ({label}) ; max_error_strain={err_strain:.3e}, max_error_stress={err_stress:.3e} -> {res}")

# 2.2. Rigid body (Translations)
# Read results from outputs generated by main.py
for i in range(3):
    axis = axes[i]
    filename = f"outputs/input_cube_trans_{axis.lower()}/input_cube_trans_{axis.lower()}_results.txt"
    
    # Simple file parser to find strain and stress
    strain_vals = []
    stress_vals = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
            
        strain_start = -1
        stress_start = -1
        for idx, line in enumerate(lines):
            if "## strain" in line:
                strain_start = idx + 1
            if "## stress" in line:
                stress_start = idx + 1
                
        if strain_start != -1:
            for line in lines[strain_start:]:
                if not line.strip(): break
                vals = [float(x) for x in line.split()]
                strain_vals.extend(vals)
                
        if stress_start != -1:
            for line in lines[stress_start:]:
                if not line.strip(): break
                vals = [float(x) for x in line.split()]
                stress_vals.extend(vals)
                
        strain_arr = np.array(strain_vals)
        stress_arr = np.array(stress_vals)
        
        if strain_arr.size == 0:
             max_strain = 999.0
             max_stress = 999.0
        else:
             max_strain = np.max(np.abs(strain_arr))
             max_stress = np.max(np.abs(stress_arr))
        
        if strain_arr.size > 0 and np.allclose(strain_arr, 0) and np.allclose(stress_arr, 0):
            res = "OK"
        else:
            res = "FAIL"
            
        print(f"MULTI ELEMENT - Rigid body (Trans {axis}) ; EXPECTED : 0.0, max_strain={max_strain:.3e} -> {res}")


# 2.3. Rigid body (Rotations)
for i in range(3):
    axis = axes[i]
    filename = f"outputs/input_cube_rot_{axis.lower()}/input_cube_rot_{axis.lower()}_results.txt"
    
    # Simple file parser to find strain and stress
    strain_vals = []
    stress_vals = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
            
        strain_start = -1
        stress_start = -1
        for idx, line in enumerate(lines):
            if "## strain" in line:
                strain_start = idx + 1
            if "## stress" in line:
                stress_start = idx + 1
                
        if strain_start != -1:
            for line in lines[strain_start:]:
                if not line.strip(): break
                vals = [float(x) for x in line.split()]
                strain_vals.extend(vals)
                
        if stress_start != -1:
            for line in lines[stress_start:]:
                if not line.strip(): break
                vals = [float(x) for x in line.split()]
                stress_vals.extend(vals)
                
        strain_arr = np.array(strain_vals)
        stress_arr = np.array(stress_vals)
        
        if strain_arr.size == 0:
             max_strain = 999.0
             max_stress = 999.0
        else:
             max_strain = np.max(np.abs(strain_arr))
             max_stress = np.max(np.abs(stress_arr))
        
        if strain_arr.size > 0 and np.allclose(strain_arr, 0) and np.allclose(stress_arr, 0):
            res = "OK"
        else:
            res = "FAIL"
            
        print(f"MULTI ELEMENT - Rigid body (Rot {axis}) ; EXPECTED : 0.0, max_strain={max_strain:.3e} -> {res}")

# 2.4. Penalty
fixed_nodes = [0,1,3,4]
loaded_nodes = [2,5,6,7]
load_val = 10.0

K_base = Ksys_cube.copy()
F_base = np.zeros(ndof_cube)
for n in loaded_nodes:
    F_base[3*n+2] = load_val

pdof_penalty = []
for n in fixed_nodes:
    for d in range(3):
        pdof_penalty.append([n+1, d+1, 0.0])
pdof_penalty = np.array(pdof_penalty)

# 2.4.1. Z too low
Z = 1.0
K_pen, F_pen = boundcond(pdof_penalty, K_base, F_base, Z, input_cube.node.shape[0])
u = np.linalg.solve(K_pen, F_pen)

max_disp_fixed = 0.0
for n in fixed_nodes:
    for d in range(3):
        val = np.abs(u[3*n+d])
        if val > max_disp_fixed:
            max_disp_fixed = val

if max_disp_fixed > 1e-3:
    res = "OK"
else:
    res = "FAIL"
print(f"MULTI ELEMENT - Penalty Z={Z:.0e} (low) ; EXPECTED : High error ; CALCULATED : {max_disp_fixed:.3f} -> {res}")

# 2.4.2 Z normal
Z = 1.0e6
K_pen, F_pen = boundcond(pdof_penalty, K_base, F_base, Z, input_cube.node.shape[0])
u_normal = np.linalg.solve(K_pen, F_pen)

max_disp_fixed = 0.0
for n in fixed_nodes:
    for d in range(3):
        val = np.abs(u_normal[3*n+d])
        if val > max_disp_fixed:
            max_disp_fixed = val

if max_disp_fixed < 1e-4:
    res = "OK"
else:
    res = "FAIL"
print(f"MULTI ELEMENT - Penalty Z={Z:.0e} (normal) ; EXPECTED : ~0.0 ; CALCULATED : {max_disp_fixed:.3e} -> {res}")

# 2.4.3. Z too high
Z = 1.0e100
K_pen, F_pen = boundcond(pdof_penalty, K_base, F_base, Z, input_cube.node.shape[0])
u_high = np.linalg.solve(K_pen, F_pen)

diff = np.max(np.abs(u_high - u_normal))

if diff > 1e-5:
    res = "OK"
else:
    res = "FAIL"
print(f"MULTI ELEMENT - Penalty Z={Z:.0e} (high) ; EXPECTED : Instability ; CALCULATED : diff_with_normal={diff:.3e} -> {res}")