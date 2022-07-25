import csv
from random import *
import re
#from Sequence_Analysis import Sequence_Length

#SAMPLE DATA BELOW

#The current sample being used is something, idk (fr this time, i have no clue)

#replace with user input (also set +/- gap for seq len)
sample_exp = {"seq_length" : 598.1,   "fragments_MseI" : [469.8, 400, 140.5],  "fragments_Hpy188I" : [598.1, 389.4, 215.8]}
sample_exp["fragments_MseI"] = [x for x in sample_exp["fragments_MseI"] if x <= sample_exp["seq_length"]]
sample_exp["fragments_Hpy188I"] = [x for x in sample_exp["fragments_Hpy188I"] if x <= sample_exp["seq_length"]]



#FORMATS SEQUENCE LIST
def db_import(db_path):
    with open(db_path, 'r') as seqFile:
        
        #creates list of definitions for all sequences
        seq = csv.reader(seqFile, delimiter='|')
        database = []
        

        #iterates through sequence list
        next(seq)
        for row in seq:

            data = {}
            
            data["defline"] = row[0]
            data["seq_length"] = int(row[1])
            data["fragments_MseI"] = [int(x) for x in row[2].split(',')]
            data["fragments_Hpy188I"] = [int(x) for x in row[3].split(',')]

            #adds the definitions to a list
            database.append(data)
            
            next(seq)
    return database

#LENGTH MATCHING
def seqlen_filter(gel_exp, gel_db, buffer=100):
    db_lenfilt = []
    for bacteria in gel_db:
        if abs(gel_exp["seq_length"]-bacteria["seq_length"]) <= buffer:
            db_lenfilt.append(bacteria)
    return db_lenfilt

#FRAGMENT NUMBER MATCHING
def fragnum_filter(gel_exp, gel_db, reject=100):
    #If more bands in experimental than digital = not a match 
    #If more bands in the digital than the experimental, may be a match
    db_fragfilt = []
    for bacteria in gel_db:
        fragments_msei = [x for x in bacteria["fragments_MseI"] if x >= reject]
        fragments_hpy188i = [x for x in bacteria["fragments_Hpy188I"] if x >= reject]
        if len(fragments_msei) >= len(gel_exp['fragments_MseI']) and len(fragments_hpy188i) >= len(gel_exp['fragments_Hpy188I']):
            db_fragfilt.append(bacteria)
    return db_fragfilt


#FRAGMENT LENGTH MATCHING (IS EXP. VALUE IN THE RANGE OF THE DB. VALUE)
def fraglen_filter_dbrange(gel_exp, gel_db, reject=100, range=100):
    db_fraglenfilt_dbrange = []

    #Iterates through all bacteria in the database
    for bacteria in gel_db:

        #removes fragments under 100 and sorts database by length
        db_msei = [x for x in bacteria["fragments_MseI"] if x >= reject]
        db_msei.sort()
        db_hpy188i = [x for x in bacteria["fragments_Hpy188I"] if x >= reject]
        db_hpy188i.sort()
        exp_msei = gel_exp["fragments_MseI"]
        exp_msei.sort()
        exp_hpy188i = gel_exp["fragments_Hpy188I"]
        exp_hpy188i.sort()
        
        #All bacteria initially assumed to be a match
        potential_match = True

        #Iterates through every database value for msei
        for dbv in db_msei:

            #Sets a range for acceptable database values
            lower_range = dbv - range
            upper_range = dbv + range

            #Experimental value initially assumed to fail the range check
            range_check = False

            #Iterates through every experimental cut value and checks to see if there is a matching database value
            for expv in exp_msei:
                if lower_range <= expv and expv <= upper_range:
                    range_check = True

            if not range_check:
                potential_match = False
                break

        #Iterates through every database value for msei
        for dbv in db_hpy188i:

            #Sets a range for acceptable database values
            lower_range = dbv - range
            upper_range = dbv + range

            #Experimental value initially assumed to fail the range check
            range_check = False

            #Iterates through every experimental cut value and checks to see if there is a matching database value
            for expv in exp_hpy188i:
                if lower_range <= expv and expv <= upper_range:
                    range_check = True
            if not range_check:
                potential_match = False
                break
                
        #If there is a potential match, it adds the database bacteria to a new list of potential matches
        if potential_match:
            db_fraglenfilt_dbrange.append(bacteria)
    return db_fraglenfilt_dbrange


def fraglen_filter_exprange(gel_exp, gel_db, range=100):
    db_fraglenfilt_exprange = []

    #Iterates through all bacteria in the database
    for bacteria in gel_db:

        #Creates lists of fragment lengths and sorts database by length
        db_msei = [x for x in bacteria["fragments_MseI"]]
        db_msei.sort()
        db_hpy188i = [x for x in bacteria["fragments_Hpy188I"]]
        db_hpy188i.sort()
        exp_msei = gel_exp["fragments_MseI"]
        exp_msei.sort()
        exp_hpy188i = gel_exp["fragments_Hpy188I"]
        exp_hpy188i.sort()
        
        db_full_length = bacteria["seq_length"]
        
        #All bacteria initially assumed to be a match
        potential_match = True

        #Iterates through every experimental value for msei
        for expv in exp_msei:

            #Sets a range for acceptable database values
            lower_range = expv - range
            upper_range = expv + range

            #Database value initially assumed to fail the range check
            range_check = False

            #Iterates through every database cut value and checks to see if there is a matching experimental value
            for dbv in db_msei:
                if lower_range <= dbv and dbv <= upper_range:
                    range_check = True


            if lower_range <= db_full_length and db_full_length <= upper_range:
                range_check = True

            #Sequence length check needed

            if not range_check:
                potential_match = False
                break

        #Iterates through every experimental value for msei
        for expv in exp_hpy188i:

            #Sets a range for acceptable database values
            lower_range = expv - range
            upper_range = expv + range

            #database value initially assumed to fail the range check
            range_check = False

            #Iterates through every database cut value and checks to see if there is a matching experimental value
            for dbv in db_hpy188i:
                if lower_range <= dbv and dbv <= upper_range:
                    range_check = True

            if lower_range <= db_full_length and db_full_length <= upper_range:
                range_check = True
            if not range_check:
                potential_match = False
                break
                
        #If there is a potential match, it adds the database bacteria to a new list of potential matches
        if potential_match:
            db_fraglenfilt_exprange.append(bacteria)
    return db_fraglenfilt_exprange



#FRAGMENT STATISTICAL ANALYSIS
# for each database length, the difference between its closest length in the sample 
# is squared and then added to the sum
def difference_value(exp_lengths, db_lengths):
    sum = 0
    # print(str(exp_lengths) + " " + str(db_lengths)) #debugging code
    # iterates through the db_lengths
    for db_length in db_lengths:
        least_difference = 99999
        # iterates through the exp_lengths
        for exp_length in exp_lengths:
            if abs(db_length-exp_length) < least_difference:
                least_difference = abs(db_length-exp_length)
        sum += least_difference**2
    return float(sum)/len(db_lengths)
    
def rank(db_fraglen_exprange, sample_exp):
    db_ranked = []
    for bacteria in db_fraglen_exprange:
        db_ranked.append(bacteria)
    for bacteria in db_ranked:
        db_diffvalHyp188I = difference_value(sample_exp["fragments_Hpy188I"], bacteria["fragments_Hpy188I"])
        db_diffvalMseI = difference_value(sample_exp["fragments_MseI"], bacteria["fragments_MseI"])
        exp_diffvalHyp188I = difference_value(bacteria["fragments_Hpy188I"], sample_exp["fragments_Hpy188I"])
        exp_diffvalMseI = difference_value(bacteria["fragments_MseI"], sample_exp["fragments_MseI"])
        db_diffval = (db_diffvalHyp188I + db_diffvalMseI + exp_diffvalHyp188I + exp_diffvalMseI)/4
        bacteria["diffval"] = db_diffval
        #print(bacteria["diffval"])
        #print(f"{i['defline']} \t {i['diffval']}")
    
    rank_sorted = sorted(db_ranked, key=lambda x: x["diffval"])
    return rank_sorted


import random

def main(passed_sample):
    db_init = db_import('./Sequence_Analyses.csv')

    random_index = random.randint(0,len(db_init) -1)
    sample_value = db_init[random_index] #change this to passed sample to use something passed in

    #print(len(db_init))
    db_lenfilt = seqlen_filter(sample_value, db_init)
    #print(len(db_lenfilt))
    #db_fragfilt = fragnum_filter(sample_exp, db_lenfilt)
    #print(len(db_fragfilt))
    db_fraglen_dbrange = fraglen_filter_dbrange(sample_value, db_lenfilt)
    #print(len(db_fraglen_dbrange))
    db_fraglen_exprange = fraglen_filter_exprange(sample_value, db_fraglen_dbrange)
    #print(len(db_fraglen_exprange))



    db_sorted_final_questionmark = rank(db_fraglen_exprange, sample_value)
    for bacteria_iter in range(0,5):
        print(db_sorted_final_questionmark[bacteria_iter]["defline"])
        print("Diffval: " + str(db_sorted_final_questionmark[bacteria_iter]["diffval"]))
        #print(db_sorted_final_questionmark[bacteria_iter])
        #break
        print("\n")
    print("Total number of matches: " + str(len(db_sorted_final_questionmark)))

    file_name = f'{random.randomint(0,10000)}_manual.txt'
    if('defline' in sample_value):
        file_name = re.match(r'S[0-9]+',sample_value["defline"])[0]
    h = open(f'./{file_name}.txt', 'w')
    for testdta in db_sorted_final_questionmark:
        h.write(sample_value["defline"] + "\n")
    h.close()

main(sample_exp)