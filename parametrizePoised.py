from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.SaltRemover import SaltRemover
import collections
from simtk.openmm.app import NoCutoff, HBonds
from simtk import openmm
from simtk import unit
import numpy as np


force_field = ForceField('test_forcefields/smirnoff99Frosst.offxml')    #load the amber compatible smirnoff FF

cntErrors1=[]   # a list to track errors when using pure openff parametrization
cntErrors2=[]   # a list to track errors when using rdkit and openff parametrization
failedMols=[]   # a collector for failed molecules

with open('data/fragment_library.smi', 'r') as reader:
    for smiles in reader.readlines():
        smiles=smiles.strip()
        try:
            mol = Molecule.from_smiles(smiles)

        except Exception as err:
            cntErrors1.append(type(err).__name__)
            pass
reader.close()


print(collections.Counter(cntErrors1))
#img=Draw.MolsToGridImage(failedMols,molsPerRow=6,subImgSize=(200,200))
#img.save('test.png')

sdfwriter = Chem.SDWriter('output/minimized.sdf')

with open('data/fragment_library.smi', 'r') as reader:
    n=0
    for smiles in reader.readlines():
        n+=1
        smiles=smiles.strip()
        print(str(n)+" : "+smiles)
        rdkit_mol=Chem.MolFromSmiles(smiles)
        rdkit_mol_h = Chem.AddHs(rdkit_mol)
        AllChem.EmbedMolecule(rdkit_mol_h)  #generate a 3D conformation to avoid the errors before
        Chem.rdmolops.AssignAtomChiralTagsFromStructure(rdkit_mol_h)    #very important, else you'll repeat the errors from before
        
        AllChem.ComputeGasteigerCharges(rdkit_mol_h)    #thats not at all mandatory nore recommended. It's just to go faster here. Don't do that in your code !!!
        
        try:
            #here is where the magic happens
            mol = Molecule.from_rdkit(rdkit_mol_h)  #create a openff molecule from the rdkit one
            top=mol.to_topology()                   #extract the topology

            ligand_system = force_field.create_openmm_system(top,charge_from_molecules=[mol])       #create our openmm system (that's the simplest version here)

            #classical setup of an openmm simulation
            #here we'll do just a bit of minimization
            integrator = openmm.LangevinIntegrator(300*unit.kelvin,91/unit.picosecond, 0.002*unit.picoseconds)
            simulation = openmm.app.Simulation(top, ligand_system, integrator)
            simulation.context.setPositions(mol.conformers[0])
            print("     starting minimization")
            state = simulation.context.getState(getEnergy=True, getForces=True)
            print('     Starting pot energy:', state.getPotentialEnergy())
            simulation.minimizeEnergy(tolerance=0, maxIterations=10000)
            state = simulation.context.getState(getPositions=True, getEnergy=True, getForces=True)
            print('     Ending pot energy:', state.getPotentialEnergy())
            newcoor = state.getPositions(asNumpy=True)
            print("     ",np.max(np.abs(newcoor - mol.conformers[0])))
            
            #now let's write back out the last conformer of the minimization for visual inspection
            mol.conformers[0]=newcoor
            rdkit_out=mol.to_rdkit()
            sdfwriter.write(rdkit_out)
            sdfwriter.flush()
            
        except Exception as err:
            print(err)
            failedMols.append(rdkit_mol_h)
            cntErrors2.append(type(err).__name__)


print(collections.Counter(cntErrors2))
