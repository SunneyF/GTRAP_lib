
from common_functions import remove_points
def write_a_txt(data,instance,types,Delta,alg,proj):
    data_unique = remove_points(data) if type(data) != tuple else data
    with open("Results\\new_output_{}_{}_{}_{}_{}_{}.txt".format(instance,Delta[0], Delta[1],Delta[2],alg,proj), "{}".format(types)) as txt_file:
        if type(data)!=tuple:
            for line in data_unique:
                txt_file.write(str(line[0]) + ',' + str(line[1]) + ',' + str(line[2]) + "\n") # works with any number of elements in a line
        else:
            txt_file.write("u= " + str(data[0]) + ',' + str(data[1]) + "\n")
