original_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]

transposed_list = list(map(list, zip(*original_list)))

print(transposed_list)
