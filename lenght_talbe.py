table = {"H":{"H":74,"F":92,"Cl":127,"Br":141,"I":161}
        ,"C":{"H":109,"C":154,"Si":186,"N":147,"O":143,"P":187,"S":181,"F":133,"Cl":177,"Br":194,"I":213}
        ,"N":{"H":101,"N":146,"P":177,"O":144,"F":139,"Cl":191,"Br":214,"I":222}
        ,"O":{"H":96,"P":160,"O":148,"S":151,"F":142,"Cl":164,"Br":172,"I":194}
        ,"Si":{"H":148,"Si":234,"S":210,"O":161,"F":156,"Cl":204,"Br":216,"I":240}
        ,"P":{"H":142,"Si":227,"P":221,"F":156,"Cl":204,"Br":222,"I":246}
        ,"S":{"H":134,"S":204,"F":158,"Cl":201,"Br":225,"I":234}
        ,"F":{"F":143,"Cl":166,"Br":178,"I":187}
        ,"Cl":{"Cl":199,"Br":214,"I":243}
        ,"Br":{"Br":228,"I":248}
        ,"I":{"I":266}
}
multi_table = {2:{"C":{"C":134,"N":127,"O":123},"N":{"N":122,"O":120},"O":{"O":121}}
              ,3:{"C":{"C":121,"N":115,"O":113},"N":{"N":110,"O":106}}}

sphere_table = {"H":['#F0F0F0',25],"C":['#000000',70],"N":['#59C4F2',65],"O":['r',60],
                "F":['g',50],"Si":['#9c9c9c',110],"P":['#FF8100',100],"S":['#FFFF00',100],
                "Cl":['#3BFF3B',100],"Br":['#831F18',115],"I":['#6229AC',140],
                "B":["#FBCEB1",90],"Be":["#8977AD",112]}
# H C N O F Si P S Cl Br I

def lenght(a,b,number):
    ratio = 0.05
    if number=="-":
        try:
            return table[a][b]*ratio
        except:pass
        try:
            return table[b][a]*ratio
        except:return 134*ratio
    else:
        number = 2 if number=="=" else 3
        try:
            return multi_table[number][a][b]*ratio
        except:pass
        try:
            return multi_table[number][b][a]*ratio
        except:return 134*ratio

def radius(atom):
    return sphere_table[atom][1] * 0.07

def color(atom):
    return sphere_table[atom][0]