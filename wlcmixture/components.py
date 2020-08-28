#Import Modules
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl
import mpltern
from mpltern.ternary.datasets import get_triangular_grid
import math
from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate


class Polymer:

    def __init__(self, A, L, la, linear_charge_density, name):
        '''
        Parameters
        ----------
        A : float
            The cross sectional area of the polymer [nm^2]
        L : float
            The length of the full polymer [nm]
        la : float
            The persistence length of the polymer [nm]
        linear_charge_density : float
            The linear charge density of the polymer [e/nm]
        name : string
            Name of the polymer for keeping track of the mixture
        '''
        self.cross_sec_area = A
        self.length = L
        self.persistence_length = la
        self.linear_charge_density = linear_charge_density
        self.name = name
        self.volume = A * L
        self.charge_density = linear_charge_density / A

    def all_properties(self):
        """
        Prints out all of the calcuable properties for the polymer
        """
        print("Polymer:", self.name)
        print(" ", "Cross Sectional Area =", self.cross_sec_area, "nm^2")
        print(" ", "Polymer Length =", self.length, "nm")
        print(" ", "Polymer Persistence Length =", self.persistence_length, "nm")
        print(" ", "Linear Charge Density =", self.linear_charge_density, "e/nm")
        print(" ", "Volume of Polymer Chain =", self.volume, "nm^3")
        print(" ", "Charge Density = ", self.charge_density, "e/nm^3")
        print(" ", "Radius of Gyration =", self.radius_of_gyration(), "nm")
        print("")

    def radius_of_gyration(self):
        """
        Calculates the radius of gyration for the polymer
        """
        kuhn_length = 2 * self.persistence_length
        N = self.length / kuhn_length

        RG_2 = (N * kuhn_length ** 2) / 6
        Rg = np.sqrt(RG_2)

        return Rg

    def structure_factor(self, k, Phi_P, SF_Type):
        """
        Calculates the structure factor of the chain using either the debye function (Gaussian) or WLC

        Parameters
        ----------
        k : float
            Wavenumber
        Phi_P : float
            Volume fraction of polymer component
        SF_Type : string
            detemines what method is used to calculate the structure factor (Gaussian or WLC)
        """
        if SF_Type == 'Gaussian':
            Rg = self.radius_of_gyration()
            x = (k * Rg) ** 2
            S_P = (Phi_P / (self.cross_sec_area * self.length)) * (2 / x ** 2) * (np.exp(-x) - 1 + x)
        if SF_Type == "WLC":
            print("This Feature is not yet added")
            S_P = 1

        return S_P

    def counter_ion(self, Ion):
        """
        Links a counter ion to the polymer creating a polymer-salt for future calculations

        Parameter
        ---------
        Ion : class
            This is a premade ion class with the opposite charge sign than the polymer
        """
        if self.charge_density * Ion.charge_density < 0:
            self.counter_ion = [Ion]
        else:
            print("Polymer and counter ion have the same charge sign")


class Ion:

    def __init__(self, charge, volume, name):
        """
        Parameters
        ----------
        charge : float
            Charge on the ion
        volume : float
            Volume of the ion
        name : string
            name of the ion
        """
        self.charge = charge
        self.volume = volume
        self.name = name
        self.charge_density = charge / volume

    def all_properties(self):
        """
        Prints out all of the calcuable properties for the ion
        """
        print("Ion:", self.name)
        print(" ", "Charge =", self.charge, "e")
        print(" ", "Volume =", self.volume, "nm^3")
        print(" ", "Charge Density =", self.charge_density, "e/nm^3")
        print("")

    def structure_factor(self, Phi_I):
        """
        Calculates the structure factor of the ion

        Parameters
        ----------
        Phi_I : float
            volume fraction of ion
        Returns
        S_I : float
            structure factor of ion
        """
        v_I = self.volume
        S_I = Phi_I / v_I  # *(2*math.pi)**3
        return S_I


class Solvent:

    def __init__(self, charge_separation, charge_valency, volume, name):
        """
        Parameters
        ----------
        charge_separation : float
            the distance between the two charges of the dipole that represents the solvent
        charge_valency : float
            charge on the dipole
        volume : float
            volume of solvent molecule
        name : string
            name of the solvent
        """
        self.charge_separation = charge_separation
        self.charge_valency = charge_valency
        self.volume = volume
        self.name = name
        self.charge_density = 0

    def all_properties(self):
        """
        Prints out all of the calcuable properties for the solvent
        """
        print("Solvent:", self.name)
        print(" ", "Charge Separation =", self.charge_separation, "nm")
        print(" ", "Charge Valency =", self.charge_valency, "e")
        print(" ", "Volume =", self.volume, "nm^3")
        print(" ", "Charge Density = ", self.charge_density, "e/nm^3")
        print("")

    def structure_factor(self, Phi_S):
        """
        Calculates the structure factor of the solvent

        Parameters
        ----------
        Phi_I : float
            volume fraction of ion
        Returns
        S_I : float
            structure factor of ion
        """
        v_S = self.volume
        M_SS = v_S ** 2
        S_S = Phi_S / v_S * M_SS  # * (2*math.pi)**3
        return S_S


class Mixture:
    def __init__(self):
        """
        Prepares the empty list to add the polymers and ions to. We can only have 1 solvent
        """
        self.polymers = []
        self.ions = []
        self.solvent = [Solvent(charge_separation=1, charge_valency=1, volume=1, name='Water')]

    def components(self):
        """
        Prints out a list of all of the components in the mixture
        """
        print('Polymers:')
        for P in self.polymers:
            print('  ', P.name)
        print('Ions:')
        for I in self.ions:
            print('  ', I.name)
        print('Solvent:')
        for S in self.solvent:
            print('  ', S.name)

    def mixture(self):
        """
        Stores a list of every component and its properties

        Returns
        -------
        mix : list
            Components stored in a list starting with polymers then ions and solvent
        """
        mix = []
        for P in self.polymers:
            mix.append(P)
        for I in self.ions:
            mix.append(I)
        for S in self.solvent:
            mix.append(S)
        return mix

    def add_polymer(self, Polymer):
        """
        Adds a polymer to the solution

        Parameter
        ---------
        Polymer : class
            Pre-built polymer class
        """
        if Polymer.type == "polymer":
            self.polymers.append(Polymer)
        else:
            print("The component you tried to add was not a polymer class")

    def add_ion(self, Ion):
        """
        Adds an ion to the solution

        Parameter
        ---------
        Ion : class
            Pre-built ion class
        """
        if Ion.type == "ion":
            self.ions.append(Ion)
        else:
            print("The component you tried to add was not an ion class")

    def add_solvent(self, Solvent):
        """
        Adds a solvent to the solution

        Parameter
        ---------
        Solvent : class
            Pre-built solvent class
        """
        if Solvent.type == 'solvent':
            self.solvent = []
            self.solvent.append(Solvent)
        else:
            print("The component you tried to add was not a solvent class")

    def input_chi(self, X):
        """
        Input the Flory-Huggins parameters as a matrix. The X[i][j] elements go in the order of polymers added
        then ions added and finally solvent.

        Parameter
        ---------
        X : array or list
            Flory-Huggins interaction parameter matrix (symmetric array)
        """
        if type(X) is list:
            X = np.array(X)

        test_symmetry = np.allclose(X, X.T, rtol=1e-05, atol=1e-08)
        if test_symmetry:
            self.chi_matrix = X
            # print ("Reminder: to add output table to display the matrix with labels")
        else:
            print("Your Flory-Huggins Interaction matrix is not symmetric. Please fix this")

    def check_charge_neutrality(self, Phi):
        """
        Checks to see if the volume fractions for the components lead to a charge neutral solution

        Parameter
        ---------
        Phi : array of floats
            Volume fractions of all of the species with the ordering of all polymers, all ions, and solvent
            in the order that they where added

        Return
        ------
        True or False
        """
        if abs(1 - np.sum(Phi)) > .000001:
            print("Volume Fractions do not sum to 1")
            return False
        else:
            total_charge = 0
            for i, C in enumerate(self.mixture()):
                component_charge = C.charge_density * Phi[i]
                total_charge = total_charge + component_charge
            if total_charge == 0:
                return True
            else:
                return False

    def gamma2(self, k, Phi, return_reduced_matrix=True, SF_Type="Gaussian", T=300):
        """
        Calculates the gamma2 matrix (n x n) or the reduced n-1 x n-1 matrix

        Parameter
        ---------
        k : float
            wavenumber
        Phi : array of floats (1 x n)
            Volume fractions of all of the species with the ordering of all polymers, all ions, and solvent
            in the order that they where added
        return_reduced_matrix : True/False
            True returns the n-1 x n-1 matrix that has been reduced using compressibility
        SF_Type : string
            Determines what function is used to calculate the structure factor (Gaussian or WLC)
        T : Float
            Temperature of the solution in Kelvin

        Return
        ------
        G : n x n array
            Gamma2 matrix
        """

        if not hasattr(self, 'chi_matrix'):
            print("You have not inputed a Flory-Huggins Interaction Matrix (self.input_chi(X))")
            n = len(self.mixture())
            return np.zeros((n, n))

        n = len(self.mixture())
        G = np.zeros((n, n))
        Q = np.zeros((n, n))
        S = np.zeros((n, n))
        Gamma2 = np.zeros((n - 1, n - 1))

        # Calculating Other Constants:
        Kb = 1.38064852E-23  # m2 kg s-2 K-1
        e0 = 8.8541878128E-12  # F/m

        b_S = self.solvent[0].charge_separation
        q_s = self.solvent[0].charge_valency
        c_S = Phi[n - 1] / self.solvent[0].volume

        epsilon = e0 + b_S * q_s * c_S / (3 * Kb * T)

        X = self.chi_matrix  # *(2*math.pi)**3

        for i, comp1 in enumerate(self.mixture()):
            Q1 = comp1.charge_density

            if comp1.type == 'polymer':
                S[i][i] = comp1.structure_factor(k, Phi[i], SF_Type)
            else:
                S[i][i] = comp1.structure_factor(Phi[i])

            for j, comp2 in enumerate(self.mixture()):
                Q2 = comp2.charge_density

                if k > 0:
                    Q[i][j] = Q1 * Q2 / (Kb * T * epsilon * k ** 2)  # *(2*math.pi)**3
                else:
                    Q[i][j] = 0
        if is_invertible(S):
            G = X + Q + LA.inv(S)
        else:
            print("Structure Factor Matrix is not invertable")

        if return_reduced_matrix:
            n = G.shape[0]
            Gamma2 = np.zeros((n - 1, n - 1))

            for i in range(n - 1):
                for j in range(n - 1):
                    Gamma2[i][j] = G[i][j] - G[i][n - 1] - G[n - 1][j] + G[n - 1][n - 1]

            return Gamma2
        else:
            return G

    def ternary(self):
        if not hasattr(self.polymers[0], 'counter_ion'):
            print("No Counter Ion")

        Phi_Matrix = generate_volume_fraction(3, .1, .001)

        # Convert Phi_Matrix to it's charge neutral composition

#
#     def reduced_gamma2(self,k,Phi,SF_Type="Gaussian"):
#         """
#         Calculates the gamma2 matrix (n-1 x n-1) since we can reduce the matrix size by one element using incompressibility

#         Parameter
#         ---------
#         k : float
#             wavenumber
#         Phi : array of floats (1 x n)
#             Volume fractions of all of the species with the ordering of all polymers, all ions, and solvent
#             in the order that they where added
#         SF_Type : string
#             Determines what function is used to calculate the structure factor (Gaussian or WLC)

#         Return
#         ------
#         G : n-1 x n-1 array
#             Gamma2 matrix
#         """
#         G = self.gamma2(k,Phi)
#         n = G.shape[0]
#         Gamma2 = np.zeros((n-1,n-1))

#         for i in range(n-1):
#             for j in range(n-1):
#                 Gamma2[i][j] = G[i][j] - G[i][n-1] - G[n-1][j] + G[n-1][n-1]

#         return Gamma2

#     def neutral_volume_fractions():
#     def ternary():
#     def density_fluctuations():
#     def check_stability():
#     def eigen_g2():




