
class Dendogram:

    def __init__( self ):
        self.root = None


    class Node:
        def __init__( self ):
            self.children = list()
            self.parent = None
            self.data = None
        
