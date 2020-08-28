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
