import csv
#SAMPLE DATA BELOW

#replace with user input (also set +/- gap for seq len)
sample_exp =  { "seq_length" : 585,   "fragments_MseI" : [511,227],  "fragments_Hpy188I" : [644,454,290]}
sample_exp["fragments_MseI"] = [x for x in sample_exp["fragments_MseI"] if x <= sample_exp["seq_length"]]
sample_exp["fragments_MseI"] = [x for x in sample_exp["fragments_MseI"] if x >= 100]

sample_exp["fragments_Hpy188I"] = [x for x in sample_exp["fragments_Hpy188I"] if x <= sample_exp["seq_length"]]
sample_exp["fragments_Hpy188I"] = [x for x in sample_exp["fragments_Hpy188I"] if x >= 100]

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
        if abs(gel_exp["seq_length"]-bacteria["seq_length"]) < buffer:
            db_lenfilt.append(bacteria)
    return db_lenfilt

#FRAGMENT NUMBER MATCHING
def fragnum_filter(gel_exp, gel_db, buffer=100):
    #If more bands in experimental than digital = not a match 
    #If more bands in the digital than the experimental, may be a match
    db_fragfilt = []
    for bacteria in gel_db:
        fragments_msei = [x for x in bacteria["fragments_MseI"] if x >= buffer]
        fragments_hpy188i = [x for x in bacteria["fragments_Hpy188I"] if x >= buffer]
        if len(fragments_msei) >= len(gel_exp['fragments_MseI']) and len(fragments_hpy188i) >= len(gel_exp['fragments_Hpy188I']):
            db_fragfilt.append(bacteria)
    return db_fragfilt


#FRAGMENT LENGTH MATCHING
def fraglen_filter(gel_exp, gel_db, buffer=100, range=50):
    db_fraglenfilt = []
    for bacteria in gel_db:
        #removes fragments under 100 and sorts database by length
        db_msei = [x for x in bacteria["fragments_MseI"] if x >= buffer]
        db_msei.sort()
        db_hpy188i = [x for x in bacteria["fragments_Hpy188I"] if x >= buffer]
        db_hpy188i.sort()
        exp_msei = gel_exp["fragments_MseI"]
        exp_hpy188i = gel_exp["fragments_Hpy188I"]


        potential_match = True
        
        for dbv in db_msei:

            lower_range = dbv - range
            upper_range = dbv + range
            range_check = False
            for expv in exp_msei:
                if lower_range < expv and expv < upper_range:
                    range_check = True
            if range_check == False:
                potential_match = False
        
        for dbv in db_hpy188i:

            lower_range = dbv - range
            upper_range = dbv + range
            range_check = False
            for expv in exp_hpy188i:
                if lower_range < expv and expv < upper_range:
                    range_check = True
            if range_check == False:
                potential_match = False
        
        
        if potential_match == True:
            db_fraglenfilt.append(bacteria)
        
    return db_fraglenfilt
    
    
db_init = db_import('./Sequence_Analyses.csv')
print(len(db_init))
db_lenfilt = seqlen_filter(sample_exp, db_init)
print(len(db_lenfilt))
db_fragfilt = fragnum_filter(sample_exp, db_lenfilt)
print(len(db_fragfilt))
db_fraglen = fraglen_filter(sample_exp, db_fragfilt)
print(len(db_fraglen))

h = open('./testData.txt', 'w')
for testdta in db_fraglen:
    h.write(testdta["defline"] + "\n")
h.close()


