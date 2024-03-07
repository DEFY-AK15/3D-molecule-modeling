import draw
from lenght_talbe import *
# H C N O F Si P S Cl Br I
#확장된 옥텟규칙 : Be B
valence_electron = {"H":1,"C":4,"N":5,"O":6,"F":7,"Si":4,"P":5,"S":6,"Cl":7, "Br":7,"I":7,
                    "B":3,"Be":2}
                    
posx = [-1,0,1,0]
posy = [0,-1,0,1]
def find_unshared_electron(lewis_dot):
    result={}
    for i in range(len(lewis_dot)):
        for j in range(len(lewis_dot[i])):
            temp = lewis_dot[i][j]
            if not temp in [" ","-","=","+"]:
                connect_electron = 0
                for k in range(4):
                    y,x = i+posy[k], j+posx[k]
                    if 0<=y<len(lewis_dot) and 0<=x<len(lewis_dot[y]):
                        temp2 = lewis_dot[y][x]
                        if temp2 == "-":
                            connect_electron+=2
                        if temp2 == "=":
                            connect_electron+=4
                        if temp2 == "+":
                            connect_electron+=6
                if temp in ["H","B","Be"]:
                    continue
                elif connect_electron<8:
                    result[(i,j)] = (8-connect_electron)//2
    return result
#retrun {(원자좌표) : 비공유 전자쌍 개수,...}


posx2 = [-2,0,2,0]
posy2 = [0,-2,0,2]
def find_midE(lewis_dot):
    result= {}
    for i in range(len(lewis_dot)):
        for j in range(len(lewis_dot[i])):
            temp = lewis_dot[i][j]
            if not temp in [" ","-","=","+"]:
                connected = []
                for k in range(4):
                    y,x = i+posy[k], j+posx[k]
                    if 0<=y<len(lewis_dot) and 0<=x<len(lewis_dot[y]):
                        if lewis_dot[y][x] in ["-","=","+"]:
                            connected.append([i+posy2[k], j+posx2[k]])
                if len(connected)>1:
                    result[(i,j)] = connected
    return result
#retrun [(중심원자 좌표): [(연결된원자좌표),..],...]


#line, tri, tetah, tri_pyra, bond
def find_structure(midE, unshared_electron):
    result = {}
    for temp in midE.keys():
        i,j = temp
        shared_E = len(midE[(i,j)])
        try:
            unshared_E = unshared_electron[(i,j)]
        except: unshared_E = 0
        if unshared_E ==0:
            if shared_E == 2:
                structure = "line"
            elif shared_E == 3:
                structure = "tri"
            elif shared_E == 4:
                structure = "tetah"
        elif unshared_E == 1:
            if shared_E == 3:
                structure = "tri_pyra"
        elif unshared_E == 2:
            if shared_E == 2:
                structure = "bond"
        result[(i,j)] = structure
    return result
#return {(중심원자 좌표) : 구조,...}

def is_connected(midE_list,pos1,pos2):
    return True if list(pos2) in midE_list[pos1] else False
def get_atom(pos,lewis_dot):
    return lewis_dot[pos[0]][pos[1]]
def get_connection(pos1,pos2,lewis_dot):
    return lewis_dot[(pos1[0]+pos2[0])//2][(pos1[1]+pos2[1])//2]


def only_twoE(lewis_dot):
    atom=[]
    for i in range(len(lewis_dot)):
        for j in range(len(lewis_dot[i])):
            temp = lewis_dot[i][j]
            if not temp in [" "]:
                if temp in ["-","=","+"]:
                    connection = temp
                else:
                    atom.append(temp)
    pos = [[0,0,0],[lenght(atom[0],atom[1],connection),0,0]]
    draw.draw(pos[0],[pos[1]])
    for i in range(2):
        draw.sphere(radius(atom[i]),pos[i],color(atom[i]))
    draw.end()


def one_midE(midE_list, midE_structure,lewis_dot):
    mainE = list(midE_list.keys())[0]
    n = len(midE_list[mainE])
    sideE_atom_list = [0 for _ in range(n)]
    sideE_point_list = [0 for _ in range(n)]
    sideE_lewis_pos = [0 for _ in range(n)]

    #sideE 원소종류
    for i in range(n):
        sideE_lewis_pos[i] = midE_list[mainE][i]
        sideE_atom_list[i] = lewis_dot[sideE_lewis_pos[i][0]][sideE_lewis_pos[i][1]]

    #중심원자 선정 & 강제 서브원자 선정
    O = [lewis_dot[mainE[0]][mainE[1]],[0,0,0],mainE]
    num = lewis_dot[(mainE[0]+sideE_lewis_pos[0][0])//2][(sideE_lewis_pos[0][1]+mainE[1])//2]
    sideE_point_list[0] = [0,0,lenght(O[0],sideE_atom_list[0],num)]

    #sideE 좌표
    if midE_structure[mainE] == "bond":
        cos,sin = -0.61235527,0.79058271
        sideE_point_list[1:] = draw.find_pos_bond(sideE_point_list[0],O[1])
    elif midE_structure[mainE] == "tri":
        cos,sin = 1,0
        sideE_point_list[1:] = draw.find_pos_tri(sideE_point_list[0],O[1])
        sideE_point_list[0],sideE_point_list[1]=sideE_point_list[1],sideE_point_list[0]
    elif midE_structure[mainE] == "tetah":
        cos,sin = 1,0
        sideE_point_list[1:] = draw.find_pos_tetah(sideE_point_list[0],O[1])
    elif midE_structure[mainE] == "tri_pyra":
        cos,sin = 1,0
        sideE_point_list = draw.find_pos_tri_pyra([0,0,2],O[1])
    elif midE_structure[mainE] == "line":
        cos,sin = 0,1
        sideE_point_list[1:] = draw.find_pos_line(sideE_point_list[0],O[1])

    for i in range(len(sideE_point_list)):
        sideE_point_list[i] = draw.rotation([0,-1,0],cos,sin,sideE_point_list[i])

    #sideE 길이보정
    for i in range(n):
        num = lewis_dot[(mainE[0]+sideE_lewis_pos[i][0])//2][(sideE_lewis_pos[i][1]+mainE[1])//2]
        sideE_point_list[i] = draw.change_distance(sideE_point_list[i],O[1],lenght(O[0],sideE_atom_list[i],num))
    
    draw.draw(O[1],sideE_point_list)

    #원자 구 생성
    atom = O[0]
    draw.sphere(radius(atom),O[1],color(atom))
    for i in range(n):
        atom = sideE_atom_list[i]
        draw.sphere(radius(atom),sideE_point_list[i],color(atom))
    
    draw.end()

def many_midE(midE_list, midE_structure,lewis_dot):
    mainEs_lewis_list = list(midE_list.keys())
    mainEs_num = len(mainEs_lewis_list)
    mainEs_point_list = [0 for _ in range(mainEs_num)]
    mainEs_atom_list = [0 for _ in range(mainEs_num)]
    toDo_list = []
    for i in range(1,mainEs_num):
        if is_connected(midE_list,mainEs_lewis_list[0],mainEs_lewis_list[i]):
            frist = i
            toDo_list.append([i,0])
            break
            
    for i in range(mainEs_num):
        for j in range(i+1,mainEs_num):
            if is_connected(midE_list,mainEs_lewis_list[i],mainEs_lewis_list[j]):
                toDo_list.append([i,j])

    for i in range(mainEs_num):
        mainEs_atom_list[i] = get_atom(mainEs_lewis_list[i],lewis_dot)

    mainEs_point_list[0] = [0,0,0]
    mainEs_point_list[frist] = [lenght(mainEs_atom_list[0],mainEs_atom_list[1],get_connection(mainEs_lewis_list[0],mainEs_lewis_list[1],lewis_dot)),0,0]
    
    for temp_list in toDo_list:
        subE_index,mainE_index = temp_list
        mainE = mainEs_lewis_list[mainE_index]
        subE = mainEs_lewis_list[subE_index]
        n = len(midE_list[mainE])
        sideE_atom_list = [0 for _ in range(n)]
        sideE_point_list = [0 for _ in range(n)]
        sideE_lewis_pos = [0 for _ in range(n)]

        temp=0
        #sideE 원소종류
        for i in range(n):
            if tuple(midE_list[mainE][i]) == mainEs_lewis_list[subE_index]:
                temp = i
            sideE_lewis_pos[i] = midE_list[mainE][i]
            sideE_atom_list[i] = lewis_dot[sideE_lewis_pos[i][0]][sideE_lewis_pos[i][1]]

        sideE_point_list[temp] = mainEs_point_list[subE_index]
        #sideE 좌표랑 루이스식의 원자랑 대응
        if midE_structure[mainE] == "bond":
            templist = draw.find_pos_bond(sideE_point_list[temp],mainEs_point_list[mainE_index])
            if temp == 0:
                sideE_point_list = [mainEs_point_list[subE_index],templist[0]]
            else:
                sideE_point_list = [templist[0],mainEs_point_list[subE_index]]


        elif midE_structure[mainE] == "tri":
            templist = draw.find_pos_tri(sideE_point_list[temp],mainEs_point_list[mainE_index])
            if temp == 0:
                sideE_point_list = [mainEs_point_list[subE_index],templist[1],templist[0]]
            elif temp == 1:
                if subE[1] - mainE[1] >0:
                    sideE_point_list  = [templist[1],mainEs_point_list[subE_index],templist[0]]
                else:
                    sideE_point_list  = [templist[0],mainEs_point_list[subE_index],templist[1]]
            elif temp == 2:
                if subE[1] - mainE[1] >0:
                    sideE_point_list  = [templist[0],templist[1],mainEs_point_list[subE_index]]
                else:
                    sideE_point_list  = [templist[1],templist[0],mainEs_point_list[subE_index]]
        
        
        elif midE_structure[mainE] == "tetah":
            templist = draw.find_pos_tetah(sideE_point_list[temp],mainEs_point_list[mainE_index])
            if temp == 0:
                sideE_point_list = [mainEs_point_list[subE_index],templist[1],templist[0],templist[2]]
            elif temp == 1:
                sideE_point_list = [templist[0],mainEs_point_list[subE_index],templist[2],templist[1]]
            elif temp == 2:
                sideE_point_list = [templist[0],templist[2],mainEs_point_list[subE_index],templist[1]]
            elif temp == 3:
                sideE_point_list = [templist[2],templist[1],templist[0],mainEs_point_list[subE_index]]
        elif midE_structure[mainE] == "tri_pyra":
            sideE_point_list = draw.find_pos_tri_pyra(sideE_point_list[temp],mainEs_point_list[mainE_index])
        
        elif midE_structure[mainE] == "line":
            templist = draw.find_pos_line(sideE_point_list[temp],mainEs_point_list[mainE_index])
            if temp == 0:
                sideE_point_list = [mainEs_point_list[subE_index],templist[0]]
            else:
                sideE_point_list = [templist[0],mainEs_point_list[subE_index]]

            
        #sideE 길이보정 //
        for i in range(n):
            num = lewis_dot[(mainE[0]+sideE_lewis_pos[i][0])//2][(sideE_lewis_pos[i][1]+mainE[1])//2]
            sideE_point_list[i] = draw.change_distance(sideE_point_list[i],mainEs_point_list[mainE_index],lenght(mainEs_atom_list[mainE_index],sideE_atom_list[i],num))

        for i in range(n):
            if i==temp: continue
            atom = sideE_atom_list[i]
            draw.sphere(radius(atom),sideE_point_list[i],color(atom))

        draw.draw(mainEs_point_list[mainE_index],sideE_point_list)

        temp_list = midE_list[mainE]
        for i in range(len(temp_list)):
            for j in range(mainEs_num):
                if tuple(temp_list[i]) ==mainEs_lewis_list[j]:
                    mainEs_point_list[j] = sideE_point_list[i]


    for i in range(mainEs_num):
        atom = mainEs_atom_list[i]
        draw.sphere(radius(atom),mainEs_point_list[i],color(atom))


    draw.end()

def run(lewis):
    midE_list = find_midE(lewis)
    if len(midE_list) == 0:
        only_twoE(lewis)
    elif len(midE_list) == 1:
        one_midE(midE_list,find_structure(midE_list,find_unshared_electron(lewis)),lewis)
    else:
        many_midE(midE_list,find_structure(midE_list,find_unshared_electron(lewis)),lewis)
