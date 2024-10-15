# defect_diffusivities

## Clone directory and create environment

```bash
git clone https://github.com/MUEXLY/defect_diffusivities
python -m venv venv
venv source/bin/activate
pip install -r requirements.txt
```

## Run LAMMPS files

First, grab any LAMMPS executable. My LAMMPS has the tabgab stuff, so soft-link to that executable:

```bash
cd lammps_data/
ln -s /home/jwjeffr/software_slurm/lammps-2Aug2023/build/lmp_gpu lmp_gpu
```

Two input files `interstitial.in` and `vacancy.in` are in `lammps_data/`. To run these for CrTaVW, replace the lines:

```
# force field
pair_style eam/alloy
pair_coeff * * FeCr.eam.alloy Fe Cr

# composition
set group all type 1
set type 1 type/ratio 2 0.07 749847
```

with:

```
# force field
mass            1 51.996
mass            2 180.9479
mass            3 50.9415
mass            4 183.84

pair_style hybrid/overlay eam/fs tabgap
pair_coeff * * eam/fs Cr-Ta-V-W.eam.fs Cr Ta V W
pair_coeff * * tabgap Cr-Ta-V-W.tabgap Cr Ta V W no yes

# composition
set type 1 type/ratio 2 0.25 123456
set type 1 type/ratio 3 0.33333333334 123456
set type 1 type/ratio 4 0.5 123456
```

and submit both slurm scripts `interstitial.slurm` and `vacancy.slurm`. The example scripts are for FeCr with 7 at. % Cr. The CrTaVW runs seem to be very slow for some reason.

## Add point defects

In the environment, there will be an executable named `add-pd`, which adds the point defect to the MD run:

```bash
add-pd lammps_data/interstitial.dump lammps_data/interstitial_with_pd.dump
add-pd lammps_data/vacancy.dump lammps_data/vacancy_with_pd.dump
```

The program can also read in and write to files zipped with `gzip`.

## Get MSD curves

Run `analysis.py` to calculate MSD curves.