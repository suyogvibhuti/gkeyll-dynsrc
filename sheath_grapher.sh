#!/bin/bash

echo "Please enter filename"
read file_name
echo "Please enter folder name"
read folder_name

# creating necessary folders
cd plot_photos
mkdir -p "$folder_name"
echo "$folder_name folder created in plot_photos folder"
cd $folder_name
mkdir -p "elc"
mkdir -p "ion"
mkdir -p "field"
cd ion
mkdir "ion_nrgd"
mkdir "ion_momd"
mkdir "ion_massd"
mkdir "ion_fvel"
mkdir "ion_all"
cd ..
cd elc
mkdir "elc_nrgd"
mkdir "elc_momd"
mkdir "elc_massd"
mkdir "elc_fvel"
mkdir "elc_all"
cd ..
cd ..
cd ..
echo "all subfolders created within $folder_name"

# activating pgkyl virtual environment
cd postgkyl
source pgkyl/bin/activate
cd ..

mv "sheath_boundary_ub.txt" "plot_photos/${folder_name}/"
mv "elc_particle_count.txt" "plot_photos/${folder_name}/"
mv "ion_particle_count.txt" "plot_photos/${folder_name}/"
mv "elc_flux_count.txt" "plot_photos/${folder_name}/"
mv "ion_flux_count.txt" "plot_photos/${folder_name}/"
mv "number_density_benchmark.txt" "plot_photos/${folder_name}/"
mv "field_potential_benchmark.txt" "plot_photos/${folder_name}/"
mv "ion_velocity_benchmark.txt" "plot_photos/${folder_name}/"
mv "elc_temp.txt" "plot_photos/${folder_name}/"

# creating all necessary plots
for i in {0..100}
do
	# ion
	op_file="ion_${i}_massd_${folder_name}.png"
	pgkyl "${file_name}-ion_${i}.gkyl" sel -c0 plot --saveas "$op_file" --dpi 200 --no_show
	mv "$op_file" "plot_photos/${folder_name}/ion/ion_massd/"

	op_file="ion_${i}_momd_${folder_name}.png"
	pgkyl "${file_name}-ion_${i}.gkyl" sel -c1 plot --saveas "$op_file" --dpi 200 --no_show
	mv "$op_file" "plot_photos/${folder_name}/ion/ion_momd/"

	op_file="ion_${i}_nrgd_${folder_name}.png"
	pgkyl "${file_name}-ion_${i}.gkyl" sel -c4 plot --saveas "$op_file" --dpi 200 --no_show
	mv "$op_file" "plot_photos/${folder_name}/ion/ion_nrgd/"

	op_file="ion_${i}_fvel_${folder_name}.png"
	pgkyl "${file_name}-ion_${i}.gkyl" sel -c1 "${file_name}-ion_${i}.gkyl" sel -c0 ev "f[0] f[1] /" plot --saveas "$op_file" --dpi 200 --no_show
	mv "$op_file" "plot_photos/${folder_name}/ion/ion_fvel/"

	op_file="ion_${i}_all_${folder_name}.png"
	pgkyl "${file_name}-ion_${i}.gkyl" plot --saveas "$op_file" --dpi 400 --no_show
	mv "$op_file" "plot_photos/${folder_name}/ion/ion_all/"

	# elc
	op_file="elc_${i}_massd_${folder_name}.png"
	pgkyl "${file_name}-elc_${i}.gkyl" sel -c0 plot --saveas "$op_file" --dpi 200 --no_show
	mv "$op_file" "plot_photos/${folder_name}/elc/elc_massd/"

	op_file="elc_${i}_momd_${folder_name}.png"
	pgkyl "${file_name}-elc_${i}.gkyl" sel -c1 plot --saveas "$op_file" --dpi 200 --no_show
	mv "$op_file" "plot_photos/${folder_name}/elc/elc_momd/"

	op_file="elc_${i}_nrgd_${folder_name}.png"
	pgkyl "${file_name}-elc_${i}.gkyl" sel -c4 plot --saveas "$op_file" --dpi 200 --no_show
	mv "$op_file" "plot_photos/${folder_name}/elc/elc_nrgd/"

	op_file="elc_${i}_fvel_${folder_name}.png"
	pgkyl "${file_name}-elc_${i}.gkyl" sel -c1 "${file_name}-elc_${i}.gkyl" sel -c0 ev "f[0] f[1] /" plot --saveas "$op_file" --dpi 200 --no_show
	mv "$op_file" "plot_photos/${folder_name}/elc/elc_fvel/"

	op_file="elc_${i}_all_${folder_name}.png"
	pgkyl "${file_name}-elc_${i}.gkyl" plot --saveas "$op_file" --dpi 400 --no_show
	mv "$op_file" "plot_photos/${folder_name}/elc/elc_all/"

	# field
	op_file="field_${i}_elec_${folder_name}.png"
	pgkyl "${file_name}-field_${i}.gkyl" sel -c0 plot --saveas "$op_file" --dpi 200 --no_show
	mv "$op_file" "plot_photos/${folder_name}/field"

	# update console
	echo "all .gkyl file plots created for save instance ${i}"
done
