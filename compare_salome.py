import numpy as np

salome_file = "Salome2Python/salome_export.csv"
python_file = "outputs/input_salome2py_spring/input_salome2py_spring_results.npz"

salome_vm = []

with open(salome_file, 'r') as f:
    lines = f.readlines()

    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) < 19:
            continue

        sxx = float(parts[13])
        syy = float(parts[14])
        szz = float(parts[15])
        sxy = float(parts[16])
        syz = float(parts[17])
        sxz = float(parts[18])

        vm = np.sqrt(0.5 * ((sxx-syy)**2 + (syy-szz)**2 + (szz-sxx)**2 + 6*(sxy**2 + syz**2 + sxz**2)))
        salome_vm.append(vm)

max_vm_salome = max(salome_vm)

data = np.load(python_file)

sigma = data['stress']

python_vm = []

for i in range(sigma.shape[0]):
    sxx = sigma[i,0]
    syy = sigma[i,1]
    szz = sigma[i,2]
    sxy = sigma[i,3]
    syz = sigma[i,4]
    sxz = sigma[i,5]

    vm = np.sqrt(0.5*((sxx-syy)**2 + (syy-szz)**2 + (szz-sxx)**2 + 6*(sxy**2 + syz**2 + sxz**2)))
    python_vm.append(vm)

max_vm_python = max(python_vm)


print(f"Salome : {max_vm_salome:.4f} MPa")
print(f"Python : {max_vm_python:.4f} MPa")

error = abs(max_vm_salome-max_vm_python) / max_vm_salome*100
print(f"Relative error : {error:.2f} %")