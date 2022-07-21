#8

def difference_value(exp_lengths, db_lengths):
    sum = 0
    if (len(exp_lengths) != len(db_lengths)):
        print("In not equal " + str(len(exp_lengths)) + " " + str(len(db_lengths)))
        for i in db_lengths:
            least_difference = 99999
            for j in exp_lengths:
                if abs(i-j) < least_difference:
                    least_difference = abs(i-j)
            sum += least_difference**2
            print(sum)
        return float(sum)/len(db_lengths)
    else:
        print("In equal")
        for i in range(len(db_lengths)):
            sum += (exp_lengths[i] - db_lengths[i])**2
            print(str(exp_lengths[i]) + " " + str(db_lengths[i]))
        return float(sum)/len(db_lengths)