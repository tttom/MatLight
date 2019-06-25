classdef Const
    properties (Constant)
        hPlanck=6.6260695729e-34; % J s, last two digits not certain
        hBarPlanck=Const.hPlanck/(2*pi); % J s
        G=6.6738480e-11; % m^3/kg s^2, last two digits not certain
        
        c=299792458; % m/s
        mu0=pi*4e-7; % N/A^2
        epsilon0=1/(Const.mu0*Const.c^2); % F/m
        
        mElectron=9.1093829140e-31; % kg, last two digits not certain
        mProton=1.67262177774e-27; % kg, last two digits not certain
        NAvogadro=6.02214078e23; % error range 0.00000018e23
        
        Rgas=8.314462175; % J /K mol, last two digits not certain
        kBoltzmann=Const.Rgas/Const.NAvogadro; % J/K
        
        waterFreezeT=273.15; % K
        
        unit=Unit(); u=Const.unit;
    end
    methods(Static)
        function tempInC=K2C(tempInK)
            tempInC=tempInK-Const.waterFreezeT;
        end
        function tempInK=C2K(tempInC)
            tempInK=tempInC+Const.waterFreezeT;
        end
        function tempInF=K2F(tempInK)
            tempInF=tempInK*9/5-459.67;
        end
        function tempInK=F2K(tempInF)
            tempInK=(tempInF+459.67)*5/9;
        end
    end
end