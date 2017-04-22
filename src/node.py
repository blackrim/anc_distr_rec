
class Node:
    def __init__(self):
        self.label = ""
        self.length = 0.0
        self.parent = None
        self.children = []
        self.data = {}
        self.istip = False
    
    def add_child(self,child):
        #make sure that the child is not already in there
        assert child not in self.children
        self.children.append(child)
        child.parent = self
    
    def remove_child(self,child):
        #make sure that the child is in there
        assert child in self.children
        self.children.remove(child)
        child.parent = None
    
    def leaves(self,v=None):
        if v == None:
            v = []
        if len(self.children) == 0:
            v.append(self)
        else:
            for child in self.children:
                child.leaves(v)
        return v

    def leaves_fancy(self):
        return [n for n in self.iternodes() if n.istip ]

    def lvsnms(self):
        return [n.label for n in self.iternodes() if n.istip ]

    def iternodes(self,order="preorder"):
        if order.lower() == "preorder":
            yield self
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order.lower() == "postorder":
            yield self
    
    def prune(self):
        p = self.parent
        if p != None:
            p.remove_child(self)
        return p

    def get_newick_repr(self,showbl=False):
        ret = ""
        for i in range(len(self.children)):
            if i == 0:
                ret += "("
            ret += self.children[i].get_newick_repr(showbl)
            if i == len(self.children)-1:
                ret += ")"
            else:
                ret += ","
        if self.label != None:
            ret += self.label
        if showbl == True:
            ret += ":" + str(self.length)
        return ret

