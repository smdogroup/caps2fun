#import statements

class NacaOMLOptimization(ParOpt):
    def __init__():
        #initialize pyCAPS problem

        #initialize each of the AIMs for the first time

    def analysis(self, x):
        #set design variables

        #generate structure mesh with egads and tacs AIMs

        #generate fluid mesh with Pointwise

        #initialize FUN3D AIM from Pointwise mesh

        #Setup TACS model for FUNtoFEM

        #Setup files to run FUN3D in FUNtoFEM

        #run FUNtoFEM forward analysis

        #run FUNtoFEM adjoint analysis

        #extract deriv w.r.t. struct DVs

        #extract deriv w.r.t. aero DVs

        #extract deriv w.r.t. aero mesh

        #extract deriv w.r.t. struct surf mesh

        #reorder the nodes after MPI

        #print the mesh sensitivities to .sens file

        #run FUN3D AIM post analysis to read .sens file

        #read FUN3D AIM dynout for functions, gradients

        #store functions, gradients

    def getFunction(self, x):

    def getGradient(self, x):


    def objCon():
        #objectives and constraints
    def objConGrad():
        #objectives and constraint gradients

#call the class and initialize it