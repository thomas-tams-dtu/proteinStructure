 

#Read cancermut file
cancermut_file = open("cancermuts_MZF1.csv", "r")

cancermut_dict = {}


for line in cancermut_file:
    if line.startswith(",U"):
        continue

    else: 
        columns = line.strip().split(",") 

        # Convert to integer and check the condition
        position = int(columns[2])  # Assuming the position is in the third column after splitting
        
        if position > 43 and position < 126:
            
            # Use the position as the key in the dictionary
            if position not in cancermut_dict:
                
                cancermut_dict[position] = [[columns[3], columns[4], columns[5]]]
    
            
            else:
                cancermut_dict[position].append([columns[3], columns[4], columns[5]])



print(cancermut_dict[74])



mutatex_energies_file = open("results/mutation_ddgs/energies.csv", "r")

mutatex_energies_dict = {}



mut_order = []

for line in mutatex_energies_file:


    #Extract order of columns - mutations
    if line.startswith("WT"):

        columns = line.strip().split(",")

        
        for i in range(len(columns)):
            if i > 2:
                mut_order.append(columns[i])

    
    #Extract energies 
    else:
        columns = line.strip().split(",")

        mutatex_energies_dict[columns[2]] = []

        
        for i in range(len(columns)):
            if i > 2:
                mutatex_energies_dict[columns[2]].append(columns[i])



print(mut_order)


for position, pairs in mutatex_energies_dict.items():
    print(f"Position {position}: {pairs}")





#Filter out all mutations with Gibbs energy above 3 and therefore causing 


threshold = 3

filtered_mutations = {}


n = 0

for position in mutatex_energies_dict.keys():

    filtered_mutations[position] = []

    n = 0

    for mutation in mutatex_energies_dict[position]:

        if float(mutation) < threshold:

            filtered_mutations[position].append([mut_order[n], mutation])
        
        n += 1



for position, pairs in filtered_mutations.items():
    print(f"Position {position}: {pairs}")
