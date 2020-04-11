counter=1

num_samples="20000"
basis_fct="squares" #trigonometric, squares, cookies or polynomial
num_basis="3 5"
mu="0.1 0.001 0.0001"
sigma="0" #only relevant for trgonometric, but always requires some value
mode="x" #only relevant for cookies, but always requires some value
part_max=5

for basis in $basis_fct; do
	for size in $num_basis; do
		for m in $mu; do
			for s in $sigma; do
				for mde in $mode; do
					if [ $basis == "trigonometric" ]; then
						folder_name="${basis}${size}_mu${m}_s${s}"
					elif [ $basis == "cookies" ]; then	
						folder_name="${basis}${size}_mu${m}_${mde}"
					else
						folder_name="${basis}${size}_mu${m}"
					fi	
					echo $folder_name
					for part in $(seq 1 $part_max); do
						echo "python gen_data.py \"{'num_samples': $num_samples, 'basis_fct': '$basis', 'part': $part,'folder_name': '$folder_name', 'mode': '$mde', 'num_basis': $size, 'mu': $m, 'sigma': $s}\"" > job_data$counter.sh
						echo "result=\$?"  >> job_data$counter.sh
						echo "if [ \$result -eq 0 ]; then"  >> job_data$counter.sh
						echo "	rm -- \"\$0\"" >> job_data$counter.sh
						echo "fi" >> job_data$counter.sh
						((counter++))
					done
				done	
			done		
		done	
	done	
done
		
