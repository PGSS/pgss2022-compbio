genus_dict = {}
genus_list = []

with open('./Experimental_Data_Matches/Lactobacillus Brevis Matches.txt') as f:
    lines = f.read().split(";")
    for pos in range(len(lines)):
        if 'genus' in lines[pos]:
            if lines[pos-1] in genus_list:
                genus_dict[lines[pos-1]] += 1
            else:
                genus_dict[lines[pos-1]] = 1
                genus_list.append(lines[pos-1])

print(genus_dict)


