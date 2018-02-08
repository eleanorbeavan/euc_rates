# without fossils
./coevol -d "alignment_file" -t "tree_file" -c "trait_file" -bd -f "output_file_name"

# with fossils
./coevol -d "alignment_file" -t "tree_file" -c "trait_file" -cal "fossil_file" "mean" "standard deviation" -bd -f "output_file_name"
