# AtomGrid Toolkit
# A robust tool for processing, analyzing, and converting 3D grid data
# from computational chemistry formats like DMol3 .grd, VASP CHGCAR, and Gaussian .cube.

import numpy as np
import os
import sys
import re
from pathlib import Path

# ==============================================================================
#  0. CONSTANTS
# ==============================================================================
BOHR_PER_ANGSTROM = 1 / 0.529177249  # Conversion factor from Angstrom to Bohr
ANGSTROM_PER_BOHR = 0.529177249      # Conversion factor from Bohr to Angstrom
ATOMIC_NUMBERS = { 'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86 }

# ==============================================================================
#  1. CORE DATA PROCESSING MODULE
# ==============================================================================

class CoreProcessor:
    """
    Encapsulates all core logic for file parsing, calculation, and writing.
    This design is robust, avoids code repetition, and ensures correctness.
    """
    @staticmethod
    def parse_grd_file(grd_filepath: Path):
        """
        Imports a DMol3 .grd file, correctly handling the N-1 grid dimension convention.
        Returns density in e/Å³ and cell parameters in Å.
        """
        if not grd_filepath.exists():
            raise FileNotFoundError(f"Data file not found: {grd_filepath}")

        print(f"  -> Reading DMol3 .grd: {grd_filepath.name}")
        with grd_filepath.open('r') as f:
            header_lines = [next(f) for _ in range(5)]

        lattice_params_str = header_lines[2].split()
        grid_dims_str = header_lines[3].split()

        def find_first_numeric_index(str_list):
            for i, item in enumerate(str_list):
                try:
                    float(item)
                    return i
                except ValueError:
                    continue
            return 0

        l_start = find_first_numeric_index(lattice_params_str)
        g_start = find_first_numeric_index(grid_dims_str)

        lattice_params = [float(x) for x in lattice_params_str[l_start:]]
        nx, ny, nz = [int(x) + 1 for x in grid_dims_str[g_start:]]

        data_1d = np.loadtxt(str(grd_filepath), skiprows=5)

        expected_points = nx * ny * nz
        if len(data_1d) < expected_points:
            raise ValueError(f"Data mismatch in {grd_filepath.name}. Expected {expected_points}, found {len(data_1d)}.")

        density_grid = data_1d[:expected_points].reshape((nx, ny, nz), order='F')
        print(f"     - Grid size detected: {nx} x {ny} x {nz}")
        return density_grid, (nx, ny, nz), lattice_params

    @staticmethod
    def parse_chgcar_file(chgcar_filepath: Path):
        """
        Parses a VASP CHGCAR file to extract the structure and physical density grid (in e/Å³).
        """
        if not chgcar_filepath.exists():
            raise FileNotFoundError(f"Data file not found: {chgcar_filepath}")

        print(f"  -> Reading VASP CHGCAR: {chgcar_filepath.name}")
        with chgcar_filepath.open('r') as f:
            lines = f.readlines()

        comment = lines[0].strip()
        scale = float(lines[1].strip())
        lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]]) * scale

        symbols = lines[5].split()
        counts = list(map(int, lines[6].split()))

        coord_start_line = 7
        if lines[coord_start_line].strip().lower().startswith(('d', 'c', 'k')):
            coord_start_line = 8

        total_atoms = sum(counts)
        direct_coords = [list(map(float, line.split())) for line in lines[coord_start_line : coord_start_line + total_atoms]]
        
        grid_data_start_line = coord_start_line + total_atoms
        while not lines[grid_data_start_line].strip():
            grid_data_start_line += 1

        nx, ny, nz = map(int, lines[grid_data_start_line].split())

        data_str = "".join(lines[grid_data_start_line + 1:])
        flat_data = np.fromstring(data_str, sep=' ')

        volume = np.linalg.det(lattice_vectors)
        if np.isclose(volume, 0):
            raise ValueError(f"Cell volume for '{chgcar_filepath.name}' is zero, cannot process file.")

        vasp_charge_data = flat_data[:nx*ny*nz].reshape((nx, ny, nz), order='F')
        density_grid = vasp_charge_data / volume # Convert from charge (e) to density (e/Å³)

        atoms_for_cube = []
        full_symbols = [s for s, c in zip(symbols, counts) for _ in range(c)]
        for symbol, d_coords in zip(full_symbols, direct_coords):
            cart_coords = np.dot(np.array(d_coords), lattice_vectors)
            atoms_for_cube.append({
                'atomic_number': ATOMIC_NUMBERS.get(symbol.capitalize(), 0),
                'charge': 0.0, 'coords_ang': cart_coords
            })

        a, b, c = [np.linalg.norm(v) for v in lattice_vectors]
        alpha = np.rad2deg(np.arccos(np.dot(lattice_vectors[1], lattice_vectors[2]) / (b * c)))
        beta  = np.rad2deg(np.arccos(np.dot(lattice_vectors[0], lattice_vectors[2]) / (a * c)))
        gamma = np.rad2deg(np.arccos(np.dot(lattice_vectors[0], lattice_vectors[1]) / (a * b)))
        cell_params = {'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta, 'gamma': gamma}

        print(f"     - Grid size detected: {nx} x {ny} x {nz}")
        return density_grid, (nx, ny, nz), cell_params, comment, symbols, counts, direct_coords, atoms_for_cube

    @staticmethod
    def parse_cube_file(cube_filepath: Path):
        """
        Parses a Gaussian .cube file. Returns density in e/Å³ and structure info in Å.
        """
        if not cube_filepath.exists():
            raise FileNotFoundError(f"Data file not found: {cube_filepath}")

        print(f"  -> Reading Gaussian .cube: {cube_filepath.name}")
        with cube_filepath.open('r') as f:
            title = f.readline().strip()
            comment = f.readline().strip()
            
            line3 = f.readline().split()
            natoms = int(line3[0])
            origin_bohr = np.array(list(map(float, line3[1:])))

            line4 = f.readline().split(); nx = int(line4[0]); vx_bohr = np.array(list(map(float, line4[1:])))
            line5 = f.readline().split(); ny = int(line5[0]); vy_bohr = np.array(list(map(float, line5[1:])))
            line6 = f.readline().split(); nz = int(line6[0]); vz_bohr = np.array(list(map(float, line6[1:])))
            
            atoms_for_cube = []
            for _ in range(natoms):
                line = f.readline().split()
                coords_bohr = np.array(list(map(float, line[2:])))
                atoms_for_cube.append({
                    'atomic_number': int(line[0]), 'charge': float(line[1]),
                    'coords_ang': coords_bohr * ANGSTROM_PER_BOHR
                })
            
            data_str = f.read()
            flat_data_bohr = np.fromstring(data_str, sep=' ')
            
        density_grid_bohr = flat_data_bohr[:nx*ny*nz].reshape((nx, ny, nz), order='C')
        density_grid_ang = density_grid_bohr * (BOHR_PER_ANGSTROM**3) # Convert e/Bohr³ to e/Å³

        # Get full lattice vectors by multiplying the voxel vector by the number of points.
        lattice_vectors_bohr = np.array([
            vx_bohr * nx,
            vy_bohr * ny,
            vz_bohr * nz
        ])
        lattice_vectors_ang = lattice_vectors_bohr * ANGSTROM_PER_BOHR

        a, b, c = [np.linalg.norm(v) for v in lattice_vectors_ang]
        # Handle potential floating point inaccuracies for near-orthogonal cells
        dot_bc = np.clip(np.dot(lattice_vectors_ang[1], lattice_vectors_ang[2]) / (b * c), -1.0, 1.0)
        dot_ac = np.clip(np.dot(lattice_vectors_ang[0], lattice_vectors_ang[2]) / (a * c), -1.0, 1.0)
        dot_ab = np.clip(np.dot(lattice_vectors_ang[0], lattice_vectors_ang[1]) / (a * b), -1.0, 1.0)
        alpha = np.rad2deg(np.arccos(dot_bc))
        beta  = np.rad2deg(np.arccos(dot_ac))
        gamma = np.rad2deg(np.arccos(dot_ab))
        cell_params = {'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta, 'gamma': gamma}

        print(f"     - Grid size detected: {nx} x {ny} x {nz}")
        return density_grid_ang, (nx, ny, nz), cell_params, title, atoms_for_cube
        
    @staticmethod
    def parse_cif_for_structure(cif_filepath: Path):
        """
        Parses a .cif file to extract structural information for CHGCAR or .cube writing.
        """
        if not cif_filepath.exists():
            raise FileNotFoundError(f"Structure file not found: {cif_filepath}")

        print(f"  -> Reading structure from .cif: {cif_filepath.name}")
        with cif_filepath.open('r') as f:
            lines = f.readlines()

        cell_params = {}; atom_data = []; in_atom_loop = False; header_map = {}
        for line in lines:
            parts = line.strip().split();
            if not parts: continue
            key = parts[0].lower()
            if key == '_cell_length_a': cell_params['a'] = float(parts[1].split('(')[0])
            elif key == '_cell_length_b': cell_params['b'] = float(parts[1].split('(')[0])
            elif key == '_cell_length_c': cell_params['c'] = float(parts[1].split('(')[0])
            elif key == '_cell_angle_alpha': cell_params['alpha'] = float(parts[1].split('(')[0])
            elif key == '_cell_angle_beta': cell_params['beta'] = float(parts[1].split('(')[0])
            elif key == '_cell_angle_gamma': cell_params['gamma'] = float(parts[1].split('(')[0])
            elif key.startswith('_atom_site_'):
                if not in_atom_loop: header_map.clear()
                in_atom_loop = True; header_map[key] = len(header_map)
            elif in_atom_loop and not line.strip().startswith(('_', 'loop_', 'data_')):
                try:
                    symbol = parts[header_map['_atom_site_type_symbol']].capitalize()
                    x = float(parts[header_map['_atom_site_fract_x']].split('(')[0])
                    y = float(parts[header_map['_atom_site_fract_y']].split('(')[0])
                    z = float(parts[header_map['_atom_site_fract_z']].split('(')[0])
                    atom_data.append({'symbol': symbol, 'coords': [x, y, z]})
                except (IndexError, KeyError, ValueError): in_atom_loop = False

        atom_groups = {s: [] for s in sorted({a['symbol'] for a in atom_data})}
        for atom in atom_data: atom_groups[atom['symbol']].append(atom['coords'])
        atom_symbols = list(atom_groups.keys()); atom_counts = [len(v) for v in atom_groups.values()]
        direct_coords = [c for s in atom_symbols for c in atom_groups[s]]
        
        lattice_vectors = CoreProcessor._get_lattice_vectors(cell_params)
        atoms_for_cube = []
        for atom in atom_data:
            fract_coords = np.array(atom['coords'])
            cart_coords = np.dot(fract_coords, lattice_vectors)
            atoms_for_cube.append({
                'atomic_number': ATOMIC_NUMBERS.get(atom['symbol'].capitalize(), 0),
                'charge': 0.0, 'coords_ang': cart_coords
            })

        print(f"     - Elements found: {list(zip(atom_symbols, atom_counts))}")
        return cif_filepath.stem, cell_params, atom_symbols, atom_counts, direct_coords, atoms_for_cube

    @staticmethod
    def _get_lattice_vectors(cell_params):
        a, b, c = cell_params['a'], cell_params['b'], cell_params['c']
        alpha, beta, gamma = np.deg2rad([cell_params['alpha'], cell_params['beta'], cell_params['gamma']])
        va = np.array([a, 0, 0]); vb = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
        vcx = c * np.cos(beta); vcy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
        vcz_sq = c**2 - vcx**2 - vcy**2
        vcz = np.sqrt(vcz_sq) if vcz_sq > 1e-9 else 0
        vc = np.array([vcx, vcy, vcz]); return np.array([va, vb, vc])

    @staticmethod
    def write_chgcar(chgcar_filepath: Path, comment, cell_params, symbols, counts, coords, density_grid):
        print(f"  -> Writing VASP CHGCAR file: {chgcar_filepath.name}")
        lattice_vectors = CoreProcessor._get_lattice_vectors(cell_params)
        volume = np.linalg.det(lattice_vectors)
        if np.isclose(volume, 0):
            print(f"\n !!! FATAL ERROR: Calculated cell volume for '{chgcar_filepath.name}' is zero. Aborting."); return
        nx, ny, nz = density_grid.shape
        with chgcar_filepath.open('w', newline='\n') as f:
            f.write(f"{comment}\n")
            f.write("   1.0\n")
            for vec in lattice_vectors: f.write(f"  {vec[0]:12.8f}{vec[1]:12.8f}{vec[2]:12.8f}\n")
            f.write("  " + "   ".join(symbols) + "\n")
            f.write("  " + "   ".join(map(str, counts)) + "\n")
            f.write("Direct\n")
            for coord in coords: f.write(f" {coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f}\n")
            f.write("\n")
            f.write(f" {nx} {ny} {nz}\n")
            vasp_charge_data = density_grid * volume # Convert density (e/Å³) to charge (e)
            flat_data = vasp_charge_data.flatten(order='F')
            for i in range(0, len(flat_data), 5):
                line_data = flat_data[i:i+5]
                f.write(" " + " ".join([f"{x:22.12E}" for x in line_data]) + "\n")
        print(f"     - Successfully created '{chgcar_filepath.name}'")

    @staticmethod
    def write_grd(grd_filepath: Path, comment, cell_params, density_grid):
        print(f"  -> Writing DMol3 .grd file: {grd_filepath.name}")
        nx, ny, nz = density_grid.shape
        grd_dims = [nx - 1, ny - 1, nz - 1]
        flat_data = density_grid.flatten(order='F')
        with grd_filepath.open('w', newline='\n') as f:
            f.write(f"{comment}\n")
            f.write("(1p,e16.9)\n")
            lp = cell_params
            f.write(f"  {lp['a']:.3f}  {lp['b']:.3f}  {lp['c']:.3f}  {lp['alpha']:.3f}  {lp['beta']:.3f}  {lp['gamma']:.3f}\n")
            f.write(f"  {grd_dims[0]}  {grd_dims[1]}  {grd_dims[2]}\n")
            f.write(f"    1    0  {grd_dims[0]}    0  {grd_dims[1]}    0  {grd_dims[2]}\n")
            np.savetxt(f, flat_data, fmt=' %.9E')
        print(f"     - Successfully created '{grd_filepath.name}'")

    @staticmethod
    def write_cube(cube_filepath: Path, title, atoms, origin_ang, grid_dims, cell_params, density_grid):
        """
        Writes a Gaussian .cube file with support for non-orthogonal cells.
        Input density is in e/Å³. All values converted to atomic units for output.
        """
        print(f"  -> Writing Gaussian .cube file: {cube_filepath.name}")
        nx, ny, nz = grid_dims
        lattice_vectors_ang = CoreProcessor._get_lattice_vectors(cell_params)
        
        # Voxel vectors are the lattice vectors divided by the number of GRID POINTS (not intervals)
        vx_ang = lattice_vectors_ang[0] / nx
        vy_ang = lattice_vectors_ang[1] / ny
        vz_ang = lattice_vectors_ang[2] / nz

        with cube_filepath.open('w', newline='\n') as f:
            f.write(f"{title}\n"); f.write("Converted by AtomGrid Toolkit (Non-Orthogonal Corrected)\n")
            
            origin_bohr = [c * BOHR_PER_ANGSTROM for c in origin_ang]
            f.write(f"{len(atoms):5d} {origin_bohr[0]:12.6f} {origin_bohr[1]:12.6f} {origin_bohr[2]:12.6f}\n")
            
            vx_bohr = [c * BOHR_PER_ANGSTROM for c in vx_ang]
            vy_bohr = [c * BOHR_PER_ANGSTROM for c in vy_ang]
            vz_bohr = [c * BOHR_PER_ANGSTROM for c in vz_ang]
            
            f.write(f"{nx:5d} {vx_bohr[0]:12.6f} {vx_bohr[1]:12.6f} {vx_bohr[2]:12.6f}\n")
            f.write(f"{ny:5d} {vy_bohr[0]:12.6f} {vy_bohr[1]:12.6f} {vy_bohr[2]:12.6f}\n")
            f.write(f"{nz:5d} {vz_bohr[0]:12.6f} {vz_bohr[1]:12.6f} {vz_bohr[2]:12.6f}\n")
            
            for atom in atoms:
                x, y, z = [c * BOHR_PER_ANGSTROM for c in atom['coords_ang']]
                f.write(f"{atom['atomic_number']:5d} {atom['charge']:12.6f} {x:12.6f} {y:12.6f} {z:12.6f}\n")
            
            density_conversion_factor = ANGSTROM_PER_BOHR**3 # Converts e/Å³ to e/Bohr³
            
            # Use C-style order for cube files (outer loop is X, inner is Z)
            flat_data_cube_order = density_grid.flatten(order='C') * density_conversion_factor
            for i in range(0, len(flat_data_cube_order), 6):
                line_data = flat_data_cube_order[i:i+6]
                f.write(" " + " ".join([f"{v: .5E}" for v in line_data]) + "\n")

        print(f"     - Successfully created '{cube_filepath.name}'")

    @staticmethod
    def calculate_1d_profiles(density_grid, cell_params, axis: str):
        """
        This function performs the core 1D profile calculations.
        It expects all inputs (density_grid, cell_params) to be in a consistent unit system.
        """
        axis = axis.lower()
        lattice_vectors = CoreProcessor._get_lattice_vectors(cell_params)
        volume = np.linalg.det(lattice_vectors)
        nx, ny, nz = density_grid.shape
        if axis == 'x': plane_axes, integration_axis, length, num_points_on_axis = (1, 2), 0, np.linalg.norm(lattice_vectors[0]), nx
        elif axis == 'y': plane_axes, integration_axis, length, num_points_on_axis = (0, 2), 1, np.linalg.norm(lattice_vectors[1]), ny
        elif axis == 'z': plane_axes, integration_axis, length, num_points_on_axis = (0, 1), 2, np.linalg.norm(lattice_vectors[2]), nz
        else: raise ValueError(f"Invalid axis '{axis}'. Currently only 'x', 'y', and 'z' are implemented.")
        
        print(f"  -> Calculating 1D profiles along the {axis.upper()}-axis...")
        area = volume / length if length > 1e-9 else 0
        delta_length = length / (num_points_on_axis - 1) if num_points_on_axis > 1 else length
        
        # This is the sum of density values on each plane.
        density_sum_on_plane = np.sum(density_grid, axis=plane_axes)
        
        # Number of grid points in the plane of integration
        num_points_in_plane = (nx * ny * nz) / num_points_on_axis
        
        # rho_avg: Plane-averaged density (e/Length^3)
        rho_avg = density_sum_on_plane / num_points_in_plane
        
        # n_slice: Total charge in a thin slice of thickness (delta_length) (e)
        voxel_volume = volume / (nx * ny * nz)
        n_slice = density_sum_on_plane * voxel_volume
        
        # LineDensity: Plane-integrated charge density (e/Length)
        line_density = n_slice / delta_length if delta_length > 1e-9 else np.zeros_like(n_slice)
        
        # n_per_S: Surface charge density (e/Area)
        n_per_S = n_slice / area if area > 1e-9 else np.zeros_like(n_slice)
        
        # ChargeDisplacement: Cumulative charge transfer (e)
        charge_displacement = np.cumsum(n_slice)
        
        axis_coords = np.linspace(0, length, num=num_points_on_axis)
        results = { "Position": axis_coords, "LineDensity": line_density, "ChargeDisplacement": charge_displacement, "rho_avg": rho_avg, "n_slice": n_slice, "n_per_S": n_per_S }
        print("     - Calculations complete.")
        return results

# ==============================================================================
#  2. USER INTERFACE AND WORKFLOWS
# ==============================================================================

class AtomGridWorkflow:
    def __init__(self):
        self.proc = CoreProcessor()

    def run_grid_calculation(self):
        print("\n--- [a] Grid Data Calculation (Add/Subtract) ---")
        print("Calculates based on input. Example: total.grd - frag1.grd + offset.grd")
        calc_input_str = input(" > Enter calculation: ").strip()
        parts = re.split('([+-])', calc_input_str.replace(' ', ''))
        if len(parts) < 3:
            print(" ! ERROR: You must provide at least two files and one operator (e.g., file1.grd - file2.grd)."); return
        base_file_path = Path(parts[0])
        result_grid, grid_dims, base_lattice_params = self.proc.parse_grd_file(base_file_path)
        for i in range(1, len(parts), 2):
            operator = parts[i]; fragment_path = Path(parts[i+1])
            frag_grid, frag_dims, _ = self.proc.parse_grd_file(fragment_path)
            if grid_dims != frag_dims:
                print(f"\n !!! FATAL ERROR: Grid dimensions mismatch between '{base_file_path.name}' {grid_dims} and '{fragment_path.name}' {frag_dims}. Aborting."); return
            if operator == '+': print(f"  -> Adding '{fragment_path.name}'..."); result_grid += frag_grid
            elif operator == '-': print(f"  -> Subtracting '{fragment_path.name}'..."); result_grid -= frag_grid
        print("\n  -> Grid calculation complete. Proceeding to output options...")
        cif_input = input(" > If atomic position information is not required, press Enter directly.\n > If it is required, please enter the .cif file path: ").strip()
        cif_provided = bool(cif_input); output_prefix = f"{base_file_path.stem}_result"
        while True:
            choice = input(" > Please select the output format:\n   1: .grd\n   2: CHGCAR\n   3: .cube\n   4: All\n   Your choice: ").strip()
            if choice in ['1', '2', '3', '4']: break
            print(" ! Invalid choice. Please enter 1, 2, 3, or 4.")
        if cif_provided:
            cif_path = Path(cif_input)
            comment, cell, symbols, counts, coords, atoms_for_cube = self.proc.parse_cif_for_structure(cif_path)
            origin_ang = [0.0, 0.0, 0.0]
            if choice in ['1', '4']: self.proc.write_grd(Path(f"{output_prefix}.grd"), f"Result for {comment}", cell, result_grid)
            if choice in ['2', '4']: self.proc.write_chgcar(Path(f"CHGCAR_{output_prefix}"), f"Result for {comment}", cell, symbols, counts, coords, result_grid)
            if choice in ['3', '4']: self.proc.write_cube(Path(f"{output_prefix}.cube"), f"Result for {comment}", atoms_for_cube, origin_ang, grid_dims, cell, result_grid)
        else:
            cell_from_grd = dict(zip(['a', 'b', 'c', 'alpha', 'beta', 'gamma'], base_lattice_params))
            comment = "Grid Calculation Result (no atomic data)"
            if choice in ['1', '4']: self.proc.write_grd(Path(f"{output_prefix}.grd"), comment, cell_from_grd, result_grid)
            if choice in ['2', '3', '4']:
                print("\n ! WARNING: CHGCAR and .cube files require atomic position information.")
                print(" ! These files cannot be generated without a .cif file. Skipping.")
                if choice not in ['1', '4']: print(" ! No files were generated.")

    def run_1d_profile_analysis(self):
        print("\n--- [b] 1D Profile Analysis ---")
        print("Accepts .grd, CHGCAR, or .cube files.")
        profile_input = input(" > Enter data file and axis (e.g., my_cdd.grd z): ").split()
        if len(profile_input) != 2: print(" ! ERROR: Invalid input format."); return
        
        data_path, axis = Path(profile_input[0]), profile_input[1]
        density_grid_ang, cell_params_ang = None, None # All parsers return data in Angstrom-based units
        
        suffix = data_path.suffix.lower()
        if suffix == '.grd':
            density_grid_ang, _, lattice_params = self.proc.parse_grd_file(data_path)
            cell_params_ang = dict(zip(['a', 'b', 'c', 'alpha', 'beta', 'gamma'], lattice_params))
        elif 'chgcar' in data_path.name.lower():
            density_grid_ang, _, cell_params_ang, _, _, _, _, _ = self.proc.parse_chgcar_file(data_path)
        elif suffix == '.cube':
            density_grid_ang, _, cell_params_ang, _, _ = self.proc.parse_cube_file(data_path)
        else:
            print(f" ! ERROR: Unrecognized file type '{data_path.suffix}'. Please use .grd, CHGCAR, or .cube."); return

        # --- Unit Selection ---
        while True:
            unit_choice = input(" > Please select the output unit system:\n   1: Angstrom (Å)\n   2: Bohr (a.u.)\n   Your choice: ").strip()
            if unit_choice in ['1', '2']: break
            print(" ! Invalid choice. Please enter 1 or 2.")
        
        if unit_choice == '1': # Angstrom
            unit_label = "A"
            density_grid_calc = density_grid_ang
            cell_params_calc = cell_params_ang
        else: # Bohr
            unit_label = "Bohr"
            # Convert density from e/Å³ to e/Bohr³
            density_grid_calc = density_grid_ang * (ANGSTROM_PER_BOHR**3)
            # Convert cell parameters from Å to Bohr
            cell_params_calc = {k: v * BOHR_PER_ANGSTROM if k in ['a','b','c'] else v for k, v in cell_params_ang.items()}

        # --- Calculation and Output ---
        results = self.proc.calculate_1d_profiles(density_grid_calc, cell_params_calc, axis)
        
        dat_filename = Path(f"{data_path.stem}_1D_profile_{axis.upper()}_{unit_label}.dat")
        print(f"  -> Writing data to '{dat_filename.name}' for plotting...")

        # Dynamic headers based on chosen unit
        header_lines = [
            f"1D Profile Analysis along {axis.upper()}-axis for {data_path.name}",
            f"Units are electrons (e) and {unit_label}.",
            f"Column 1: Position ({unit_label})",
            f"Column 2: LineDensity (e/{unit_label})",
            f"Column 3: ChargeDisplacement (e)",
            f"Column 4: rho_avg (e/{unit_label}^3)",
            f"Column 5: n_slice (e)",
            f"Column 6: n_per_S (e/{unit_label}^2)"
        ]
        
        col_keys = list(results.keys())
        col_units = [f"({unit_label})", f"(e/{unit_label})", "(e)", f"(e/{unit_label}^3)", "(e)", f"(e/{unit_label}^2)"]
        col_headers = "".join([f"{key}{unit:<{20-len(key)}}" for key, unit in zip(col_keys, col_units)])
        
        header = "\n".join(header_lines) + f"\n{col_headers}"
        data_to_save = np.vstack(list(results.values())).T
        np.savetxt(str(dat_filename), data_to_save, fmt='%18.8f', header=header, comments='# ')
        print(f"     - Done. You can now plot Col 1 vs Col 2 for the plane-integrated profile and Col 1 vs Col 3 for the charge displacement curve.")

    def run_universal_conversion(self):
        print("\n--- [c] Universal Format Conversion (.grd, CHGCAR, .cube) ---")
        filepath = Path(input(" > Enter the file to convert: "))
        if not filepath.exists(): raise FileNotFoundError(f"File not found: {filepath}")

        # --- 1. Parse Input File ---
        density_grid, grid_dims, cell_params, title = None, None, None, filepath.stem
        atoms, symbols, counts, direct_coords = [], [], [], []
        struct_info_available = False
        
        suffix = filepath.suffix.lower()
        if suffix == '.grd':
            density_grid, grid_dims, grd_lp = self.proc.parse_grd_file(filepath)
            cell_params = dict(zip(['a', 'b', 'c', 'alpha', 'beta', 'gamma'], grd_lp))
        elif 'chgcar' in filepath.name.lower():
            density_grid, grid_dims, cell_params, title, symbols, counts, direct_coords, atoms = self.proc.parse_chgcar_file(filepath)
            struct_info_available = True
        elif suffix == '.cube':
            density_grid, grid_dims, cell_params, title, atoms = self.proc.parse_cube_file(filepath)
            struct_info_available = True
        else:
            print(f" ! ERROR: Unrecognized input file type '{suffix}'."); return

        # --- 2. Ask for Output Format ---
        while True:
            out_choice = input(" > Please select the output format:\n   1: CHGCAR\n   2: .cube\n   3: .grd\n   4: All\n   Your choice: ").strip()
            if out_choice in ['1', '2', '3', '4']: break
            print(" ! Invalid choice. Please enter 1, 2, 3, or 4.")
        
        formats_to_write = []
        if out_choice == '1': formats_to_write = ['chgcar']
        elif out_choice == '2': formats_to_write = ['cube']
        elif out_choice == '3': formats_to_write = ['grd']
        elif out_choice == '4': formats_to_write = ['chgcar', 'cube', 'grd']

        # --- 3. Ask for Optional CIF file ---
        cif_input = input(" > If you don't need to add/overwrite atomic info, press Enter.\n > To add/overwrite atomic info, enter the .cif file path: ").strip()
        if cif_input:
            cif_path = Path(cif_input)
            title, cell_params, symbols, counts, direct_coords, atoms = self.proc.parse_cif_for_structure(cif_path)
            struct_info_available = True
            print("     - Structural information from .cif will be used for all outputs.")

        # --- 4. Write Output Files ---
        output_prefix = filepath.stem
        files_written = 0
        
        if 'chgcar' in formats_to_write:
            if not struct_info_available:
                print("\n ! WARNING: CHGCAR requires atomic information. Skipping CHGCAR generation.")
            else:
                self.proc.write_chgcar(Path(f"CHGCAR_conv_{output_prefix}"), title, cell_params, symbols, counts, direct_coords, density_grid)
                files_written += 1
        
        if 'cube' in formats_to_write:
            if not struct_info_available:
                print("\n ! WARNING: .cube format requires atomic information. Skipping .cube generation.")
            else:
                origin_ang = [0.0, 0.0, 0.0] # Standard assumption for conversions
                self.proc.write_cube(Path(f"{output_prefix}_conv.cube"), title, atoms, origin_ang, grid_dims, cell_params, density_grid)
                files_written += 1

        if 'grd' in formats_to_write:
            self.proc.write_grd(Path(f"{output_prefix}_conv.grd"), title, cell_params, density_grid)
            files_written += 1
            
        if not files_written:
             print("\n ! No files were generated. This may be due to missing structural information for CHGCAR/.cube formats.")

# ==============================================================================
#  3. EXPLANATIONS AND MAIN EXECUTION
# ==============================================================================

def explain_1d_profile_quantities():
    # This function remains unchanged.
    print("\n" + "-"*80); print("Explanation of 1D Profile Quantities (Function [b])".center(80)); print("-" * 80)
    print("These quantities are standard for analyzing charge transfer in materials science.")
    print("\n[Primary Curves for Publication]")
    print("  1. Line Density (e/Å), aka Δρ(z): The plane-integrated charge density. It shows charge accumulation/depletion along an axis.");
    print("  2. Charge Displacement (e), aka ΔQ(z): The net charge transferred across a plane at position z. Quantifies total charge movement.")
    print("\n[Intermediate & Other Quantities]")
    print("  3. Plane-Averaged Density (e/Å³), aka ρ_avg(z): The average volumetric density over the XY plane.")
    print("  4. Charge in Slice (e): Total charge in a thin slice (LineDensity * Δz).")
    print("  5. Surface Charge Density (e/Å²): Charge per unit area at a given slice (n_slice / Area_plane).")
    print("-" * 80)

def display_main_menu():
    print("\n" + "="*50); print(" AtomGrid Toolkit ".center(50, "=")); print("="*50)
    print("Please choose a function:")
    print("  [a] Grid Data Calculation (Add/Subtract)")
    print("  [b] 1D Profile Analysis (for Publication Plots)")
    print("  [c] Universal Format Conversion (.grd, CHGCAR, .cube)")
    print("\n--- Help & Info ---")
    print("  [h] Help: Explain 1D Profile quantities")
    print("  [q] Quit")
    print("-"*50)
    return input(" > Your choice: ").lower().strip()

if __name__ == '__main__':
    workflow = AtomGridWorkflow()
    while True:
        choice = display_main_menu()
        try:
            if choice == 'a': workflow.run_grid_calculation()
            elif choice == 'b': workflow.run_1d_profile_analysis()
            elif choice == 'c': workflow.run_universal_conversion()
            elif choice == 'h': explain_1d_profile_quantities()
            elif choice == 'q': print("Exiting AtomGrid Toolkit. Goodbye!"); break
            else: print(" ! Invalid choice. Please try again.")
        except FileNotFoundError as e:
            print(f"\n !!! FILE ERROR: {e}. Please check the file path and ensure the file exists.")
        except Exception as e:
            print(f"\n !!! AN UNEXPECTED ERROR OCCURRED !!!")
            print(f" ! ERROR TYPE: {type(e).__name__}")
            print(f" ! MESSAGE: {e}")
            print(" ! Please check your input files and try again.")

