#!/usr/bin/env python


r"""

====================
Example Data factory
====================

.. moduleauthor:: Antonia Mey <antonia.mey@fu-berlin.de>

"""


#imports
import numpy as np
from scipy.integrate import quad
import math
import os



#some helper functions/classes
class AssymetricDoubleWellPotential ( object ):
    r"""class for an assymetic double well potential
    """
    def __init__( self ):
        self.Z = []
        self.inner_edges = np.linspace(-0.4,4.2,20)
        self.n_bins = self.inner_edges.shape[0]+1
        self.bin_centers = np.zeros(self.n_bins)
        self.bin_centers[1:-1] = self.inner_edges[1:]-0.5*( self.inner_edges[1]-self.inner_edges[0] )
        self.bin_centers[0] = -0.7
        self.bin_centers[-1] = 4.5

    def energy( self, x ):
        return 2*(x-2)-6*(x-2)**2 + (x-2)**4
    def gradient( self, x ):
        return 4*x**3-24*x**2+36*x-6
    def integ( self , x , kT ):
        return np.exp( - self.energy( x )/kT )
    def get_partition_function ( self, kT ):
        for T in kT:
            self.Z.append( quad(self.integ, -100.0, 100.0, args=(T,))[0] )
        self.Z = np.array( self.Z )
        return self.Z

class HarmonicRestraint( object ):
    r"""class for harmonic restraints used in US
    """

    def __init__( self, r0, k ):
        self.r0 = r0
        self.k = k
    def energy( self, r):
        return 0.5 * self.k * ( r - self.r0 )**2
    def gradient( self, r ):
        return self.k*(r-self.r0)


class BrownianIntegrator( object ):
    r"""
    class that does brownian dynamics integration
    """
    def __init__( self, potential, dt, beta=1.0, mass=1.0, damping=1.0 ):
        r""" constructor function
        
        Parameters
        ----------
        potential : FoldingPotential object
            object that contains all the information about the potential in which the particle svolves
        dt : double
            timestep of the integrator
        beta : double
            inverse temperature of the simulation Default = 1
        mass : double
            mass of the particle, Default = 1.0  
        damping : double
            damping constant, Default = 1.0
        """
        self.potential = potential
        # store parameters
        self.dt = dt
        self.beta = beta
        self.mass = mass
        self.damping = damping
        # compute coefficients
        self.coeff_A = dt / ( mass * damping )
        self.coeff_B = np.sqrt( 2.0 * dt / ( beta * mass * damping ) )
        self.x = None
        self.t_index = None #thermodynamic index
    def step( self, restraint = None ):
        r"""function that carries out a single integration step
        """
        gradient = self.potential.gradient( self.x )
        if None != restraint:
            gradient += restraint.gradient( self.x )
        self.x = self.x - self.coeff_A * gradient + self.coeff_B * np.random.normal( size=self.n_dim )
        pos = self.x[0]
        if self.n_dim > 1 :
            pos = np.linalg.norm(self.x)
        return np.array( [ pos , self.t_index, self.potential.energy( pos ) ] )

    def set_temperature( self, beta ):
        r"""function that sets a new inverse temperature
        Parameters
        ----------
        beta : double
            inverse temperature beta
        """
        self.beta = beta
        self.coeff_B = np.sqrt( 2.0 * self.dt / ( self.beta * self.mass * self.damping ) )
    def set_position( self, x ):
        r"""function that sets the current position of the integrator
        Parameters
        ----------
        x : double
            current position of the particle
        """
        self.x = x
        self.n_dim = x.shape[0]
    def set_t_index( self, t ):
        r"""function that sets the temperature index, should there be more than one
        Parameters
        ----------
        t : int
            temperature index in the temperature array
        """
        self.t_index = t

class STReplica( object ):
    def __init__( self,Z, integrator, kT ):
        r""" Constructor function
        Parameters
        ----------

        Z : double array
            containing normalisation factors used for weighting different temperatures
        integrator : BrownianIntegrator object
            the actual simulation
        kT : double array
            array contianing the possible temperatures at which sampling occurs
        """
        self.Z = Z
        self.integrator = integrator
        self.kT = kT
        self.trajectory = []
    def run( self, nsteps=1 ):
        r""" function that runs the replica exchange simulation
        Parameters
        ----------
        nsteps : int
            number of steps the integrator should compute between exchanges
            Default=1
        """
        for i in range( nsteps ):
            self.trajectory.append(self.integrator.step())

    def change_temperature( self ):
        r"""Metropolis hastings steps that allows to change between two randomly chosen temperatures
        """
        r = np.random.randint( self.kT.shape[0] )
        beta_new = 1.0/self.kT[r]
        beta_old = self.integrator.beta
        beta_int = np.where( self.kT == 1.0/beta_old )
        deltaG = -np.log ( self.integrator.potential.Z[r] ) + np.log( self.integrator.potential.Z[beta_int] )
        enExp = -self.trajectory[-1][2]*(beta_new-beta_old)
        exponent = enExp + deltaG
        
        if ( exponent >= 0 ) or ( np.random.random()< np.exp(exponent) ):
            self.integrator.set_temperature(1.0/self.kT[r])
            self.integrator.set_t_index(r)
            #print "New temperature is: "+str(self.integrator.beta)


class USReplica( object ):
    r""" class descriptor
    """
    def __init__( self,integrator, RS ):
        self.integrator = integrator
        self.RS = RS
        self.trajectory = []
        self.n_therm_states = len(RS)
        self.ndim=1

    def run( self, nsteps ):
        count = 0
        for r in self.RS:
            r_traj = []
            self.integrator.set_position(np.ones(self.ndim)*r.r0/np.sqrt(self.ndim))
            self.integrator.t_index = count
            for i in range( nsteps ):
                int_return = self.integrator.step(r)
                bias = self.calculate_bias(self.integrator.x)
                r_traj.append(np.hstack((int_return,bias)))
            count = count+1
            self.trajectory.append(r_traj)
        print ('US simulation completed sucessfully!')

    def calculate_bias( self, pos ):
        bias = np.zeros(self.n_therm_states)
        for b in range(self.n_therm_states):
            bias[b]=self.RS[b].energy(pos)
        return bias


def discretize( x, inner_edges ):
    if x<inner_edges[0]:
        return 0
    if x>=inner_edges[-1]:
        return inner_edges.shape[0]
    for i in range( inner_edges.shape[0]-1 ):
        if ( inner_edges[i]<=x) and ( x< inner_edges[i+1] ):
            return i+1

#the core run functions




def run_st_simulation():
    N_EXCHANGES=2000
    #checking if directories for writing exist
    directory="ST/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    kT = np.array([ 2.0, 4.0, 7.0, 10.0, 15.0])
    n_therm_states = len(kT)
    initial_x = np.array([3.5])
    initial_t = 0

    dwp = AssymetricDoubleWellPotential()
    Z = dwp.get_partition_function(kT)
    integrator = BrownianIntegrator( dwp, 0.005, 1.0/kT[initial_t], 1.0, 1.0 )
    integrator.set_t_index(initial_t)
    integrator.set_position(initial_x)
    replica = STReplica(Z, integrator, kT )
    for i in range(N_EXCHANGES):
       replica.run(100)
       replica.change_temperature()
    traj = np.array(replica.trajectory)

    n_traj_frames =  np.shape(replica.trajectory)[0] 
    fh = open( directory+"Traj.dat", 'w' )
    for t in range( n_traj_frames ):
       fh.write( "%6d %6d %+.6e" % ( discretize( traj[t,0], dwp.inner_edges ), traj[t,1], traj[t,2] / kT[traj[t,1]] ) )
       fh.write( "\n" )
    fh.close()
    np.savetxt(directory+"kT.dat", kT)

    #constructing wham file
    # -> themodynamic state
    # |
    # v markov state
    #we need a target temperature \beta_0
    target = 0
    wham_f = os.path.join(directory,"b_K_i.dat")
    fh = open(wham_f, 'w')
    for c in range( dwp.bin_centers.shape[0] ):
       for i in range( n_therm_states ):
           fh.write( " %+.6e" % ( dwp.energy( dwp.bin_centers[c] ) * (1.0/kT[i]-1.0/kT[target] ) ) )
       fh.write( "\n" )
    fh.close()

    #exact probability distribution for comparison
    fh = open( directory+"exact.dat", 'w' )
    for c in range( dwp.bin_centers.shape[0] ):
       p=np.exp(-dwp.energy( dwp.bin_centers[c] )* 1.0/kT[target])
       fh.write( "%4f %+.6e" % ( dwp.bin_centers[c] , p ) )
       fh.write("\n")
    fh.close()



def run_us_simulation():
    #checking if directories for writing exist
    directory="US/"
    if not os.path.exists(directory):
       os.makedirs(directory)

    #setting the simulation temperature
    kT = np.array( [1.0] )
    nsteps=1000

    dwp = AssymetricDoubleWellPotential()
    Z = dwp.get_partition_function(kT)
    #exact probability distribution for comparison
    e_file = os.path.join(directory,"exact.dat")
    fh = open( e_file, 'w' )
    for c in range( dwp.bin_centers.shape[0] ):
       p=np.exp(-dwp.energy( dwp.bin_centers[c] )* 1.0/kT[0])
       fh.write( "%4f %+.6e" % ( dwp.bin_centers[c] , p ) )
       fh.write("\n")
    fh.close()

    integrator = BrownianIntegrator( dwp, 0.005, 1.0/kT[0], 1.0, 1.0) 
    restraints_pos = np.linspace(-0.9,4.7,30)
    restraint_k = np.ones(30)*90
    #restraints_pos = np.linspace(-1.8,1.8,16)
    n_therm_states = restraints_pos.shape[0]


    restraints = []
    for i in range(restraints_pos.shape[0]):
       restraints.append(HarmonicRestraint(restraints_pos[i],restraint_k[i]))
    replica = USReplica(integrator, restraints)
    replica.run(nsteps)
    for r in range(restraints_pos.shape[0]):
       traj = np.array(replica.trajectory[r])
       n_traj_frames = traj.shape[0]
       r_file = os.path.join(directory,"Traj"+str(r)+".dat")
       fh =open(r_file, 'w')
       for t in range( n_traj_frames ):
           fh.write( "%6d %6d " % (discretize( traj[t,0], dwp.inner_edges ) , traj[t,1] ) )
           for j in range(restraints_pos.shape[0]):
               fh.write("%+.6e " % (traj[t,3+j]/kT[0]))
           fh.write("\n")
       fh.close()

    wham_f = os.path.join(directory,"b_K_i.dat")
    fh = open(wham_f, 'w')
    for c in range( dwp.bin_centers.shape[0] ):
       for i in range( n_therm_states ):
           fh.write( " %+.8e" % ( restraints[i].energy( dwp.bin_centers[c] ) / kT[0] ) )
       fh.write( "\n" )
    fh.close()

