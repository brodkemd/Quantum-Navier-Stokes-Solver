import os

def lines_of_file_to_matrix(lines):
    return_matrix = []
    
    for item in lines:
        # removes \n from end, then whitespace
        item = item[:-1].strip().split(' ')

        # removes empty elements
        while item.count('') != 0:
            item.remove('')
        
        new_item = []
        for num in item:
            # string to float
            new_item.append(float(num))
    

        return_matrix.append(new_item)

    return return_matrix

