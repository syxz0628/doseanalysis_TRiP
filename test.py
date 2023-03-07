def fun_append_listonebyone(write2list, datalist):
    for list2write in datalist:
        write2list.append(list2write)
    return write2list

list=['1','2']
list2=['5','6']
#[list2.append(i) for i in list]
list2=fun_append_listonebyone(list2,list)
print(list2)


