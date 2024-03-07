

def merge_lists(list1, list2):
    merged_list = list1 + list2
    return merged_list


def func2(lista, listb):
    newlist = []
    listone = 0
    listtwo = 0
    while listone < len(lista) and listtwo < len(listb):
        if lista[listone] < listb[listtwo]:
            newlist.append(lista[listone])
            listone += 1
        else:
            newlist.append(listb[listtwo])
            listtwo += 1
    newlist.extend(lista[listone:])
    newlist.extend(listb[listtwo:])

    return newlist

a = [1, 3, 5]
b = [2, 4, 6]

answer1 = merge_lists(a, b)
answer2 = func2(a, b)

print(answer1)
print(answer2)
