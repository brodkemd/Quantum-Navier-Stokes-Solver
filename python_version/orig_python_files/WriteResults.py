def WriteResults(n, Tot_X_Pts, d, U2, Mach_D, Mach_E, Mrho_D, Mrho_E, Press_D, Press_E, Temp_D, Temp_E, Rel_MachErr, Rel_MrhoErr, Rel_PressErr, Rel_TempErr, AvRelTempErr, AvPlusSDevRelTempErr, AvMinusSDevRelTempErr,AvRelMachErr, AvPlusSDevRelMachErr, AvMinusSDevRelMachErr,AvRelMrhoErr, AvPlusSDevRelMrhoErr, AvMinusSDevRelMrhoErr,AvRelPressErr, AvPlusSDevRelPressErr, AvMinusSDevRelPressErr,AvU2, ff0_throat, ff1_throat, ff2_throat):
    #WRITERESULTS writes results of Navier-Stokes solution to files.
    #   WriteResults writes primary flow variables and related results to files
    #       for eventual plotting.
    #
    #   INPUTS:
    #       n = number of subintervals used in time partition
    #       Tot_X_Pts = total number of spatial grid-points used
    #       d = number of components of ODE solution z(t)
    #       U2 = Tot_X_Pts x (n+1) array that stores calculated mass flow rate
    #               at beginning of each time subinterval, and the mass flow
    #               rate at the final time of the ODE integration.
    #       Mach_D = 1 x Tot_X_Pts array storing calculated Mach number at each
    #                   grid-point at end of ODE integration
    #       Mach_E = 1 x Tot_X_Pts array storing exact steady-state Mach 
    #                   number at each grid-point at end of ODE integration
    #       Mrho_D = 1 x Tot_X_Pts array storing calculated (dimensionless)
    #                   mass density at each grid-point at end of ODE 
    #                   integration
    #       Mrho_E = 1 x Tot_X_Pts array storing exact (dimensionless) mass 
    #                   density at each grid-point at end of ODE integration
    #       Press_D = 1 x Tot_X_Pts array storing calculated (dimensionless)
    #                   pressure at each grid-point at end of ODE integration
    #       Press_E = 1 x Tot_X_Pts array storing exact (dimensionless)
    #                   pressure at each grid-point at end of ODE integration
    #       Temp_D = 1 x Tot_X_Pts array storing calculated (dimensionless)
    #                   temperature at each grid-point at end of ODE 
    #                   integration
    #       Temp_E = 1 x Tot_X_Pts array storing exact (dimensionless)
    #                   temperature at each grid-point at end of ODE 
    #                   integration
    #       Rel_MachErr = 1 x Tot_X_Pts array storing relative error in Mach
    #                       number at end of ODE integration
    #       Rel_MrhoErr = 1 x Tot_X_Pts array storing relative error in 
    #                       dimensionless mass density at end of ODE 
    #                       integration
    #       Rel_PressErr = 1 x Tot_X_Pts array storing relative error in 
    #                       dimensionless pressure at end of ODE 
    #                       integration
    #       Rel_TempErr = 1 x Tot_X_Pts array storing relative error in 
    #                       dimensionless temperature at end of ODE 
    #                       integration
    #       AvRelTempErr = 1 x Tot_X_Pts array storing average relative
    #                       temperature error, same value to each grid-point
    #       AvPlusSDevRelTempErr = 1 x Tot_X_Pts array storing (average plus
    #                       standard deviation) relative temperature error, 
    #                       same value to each grid-point
    #       AvMinusSDevRelTempErr = 1 x Tot_X_Pts array storing (average minus
    #                       standard deviation) relative temperature error, 
    #                       same value to each grid-point
    #       AvRelMachErr = 1 x Tot_X_Pts array storing average relative Mach
    #                       number error, same value to each grid-point
    #       AvPlusSDevRelMachErr = 1 x Tot_X_Pts array storing (average plus
    #                       standard deviation) relative Mach number error, 
    #                       same value to each grid-point
    #       AvMinusSDevRelMachErr = 1 x Tot_X_Pts array storing (average minus
    #                       standard deviation) relative Mach number error, 
    #                       same value to each grid-point
    #       AvRelMrhoErr = 1 x Tot_X_Pts array storing average relative mass
    #                       density error, same value to each grid-point
    #       AvPlusSDevRelMrhoErr = 1 x Tot_X_Pts array storing (average plus
    #                       standard deviation) relative mass density error, 
    #                       same value to each grid-point
    #       AvMinusSDevRelMrhoErr = 1 x Tot_X_Pts array storing (average minus
    #                       standard deviation) relative mass density error, 
    #                       same value to each grid-point
    #       AvRelPressErr = 1 x Tot_X_Pts array storing average relative 
    #                       pressure error, same value to each grid-point
    #       AvPlusSDevRelPressErr = 1 x Tot_X_Pts array storing (average plus
    #                       standard deviation) relative pressure error, 
    #                       same value to each grid-point
    #       AvMinusSDevRelPressErr = 1 x Tot_X_Pts array storing (average minus
    #                       standard deviation) relative pressure error, 
    #                       same value to each grid-point
    #       AvU2 = 1 x Tot_X_Pts array storing average mass flow rate (over 
    #                 nozzle) at final time of ODE integration at each grid-pt
    #       ff0_throat = d x n array storing driver function value at throat at
    #                       end of each subinterval
    #       ff1_throat = d x n array storing first derivative of driver 
    #                       function at throat at end of each subinterval
    #       ff2_throat = d x n array storing second derivative of driver 
    #                       function at throat at end of each subinterval

    # write results for various flow variables to file:

    '''
    filenameU2 = open('U2vals','w')
    filenameAvU2 = open('AvU2vals','w')
    filenameMachD = open('MachDvals', 'w')
    filenameMachE = open('MachEvals','w')
    filenameMrhoD = open('MrhoDvals','w')
    filenameMrhoE = open('Mrho_Evals','w')
    filenamePressD = open('PressDvals','w')
    filenamePressE = open('PressEvals','w')
    filenameTempD = open('TempDvals','w')
    filenameTempE = open('TempEvals','w')
    
    column = n+1
    row = Tot_X_Pts

    for gridpt in range(row + 1):
        fprintf(filenameAvU2, '#8.3f', AvU2(gridpt))
        
        if gridpt == 1
            for tyme = 1:column
                fprintf(filenameU2, '#8.3f', U2(gridpt, tyme))
            end
        elseif gridpt ~= 1
            fprintf(filenameU2, '\n')
            
            for tyme = 1:column
                fprintf(filenameU2, '#8.3f', U2(gridpt, tyme))
            end
        end
    end

    for gridpt = 1:Tot_X_Pts
        fprintf(filenameMachD, '#6.3f', Mach_D(gridpt))
        fprintf(filenameMachE, '#6.3f', Mach_E(gridpt))
        fprintf(filenameMrhoD, '#6.3f', Mrho_D(gridpt))
        fprintf(filenameMrhoE, '#6.3f', Mrho_E(gridpt))
        fprintf(filenamePressD, '#6.3f', Press_D(gridpt))
        fprintf(filenamePressE, '#6.3f', Press_E(gridpt))
        fprintf(filenameTempD, '#6.3f', Temp_D(gridpt))
        fprintf(filenameTempE, '#6.3f', Temp_E(gridpt))
    end

    filenameU2)
    filenameAvU2)
    (filenameMachD)
    (filenameMachE)
    (filenameMrhoD)
    (filenameMrhoE)
    (filenamePressD)
    (filenamePressE)
    (filenameTempD)
    (filenameTempE)

    # open files for writing flow variable relative errors

    filenameRelMachErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\RelMachErrvals',...
    'w+')
    filenameRelMrhoErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\RelMrhoErrvals',...
    'w+')
    filenameRelPressErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\RelPressErrvals',...
    'w+')
    filenameRelTempErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\RelTempErrvals',...
    'w+')

    # write data to files

    for gridpt = 1:Tot_X_Pts
        fprintf(filenameRelMachErr, '#6.3f', Rel_MachErr(gridpt))
        fprintf(filenameRelMrhoErr, '#6.3f', Rel_MrhoErr(gridpt))
        fprintf(filenameRelPressErr, '#6.3f', Rel_PressErr(gridpt))
        fprintf(filenameRelTempErr, '#6.3f', Rel_TempErr(gridpt))
    end

    # close files 

    (filenameRelMachErr)
    (filenameRelMrhoErr)
    (filenameRelPressErr)
    (filenameRelTempErr)

    # open files to write mean and mean +/- standard deviation of relative
    #   errors

    filenameAvRelTempErr = open(...
        'C:\Users\Marek\OneDrive\Documents\ns_q\AvRelTempErr', ...
        'w+')
    filenameAvRelMachErr = open(...
        'C:\Users\Marek\OneDrive\Documents\ns_q\AvRelMachErr', ...
        'w+')
    filenameAvRelMrhoErr = open(...
        'C:\Users\Marek\OneDrive\Documents\ns_q\AvRelMrhoErr', ...
        'w+')
    filenameAvRelPressErr = open(...
        'C:\Users\Marek\OneDrive\Documents\ns_q\AvRelPressErr', ...
        'w+')
    
    filenameAvPlusSDevRelTempErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\AvPlusSDevRelTempErr',...
        'w+')
    filenameAvPlusSDevRelMachErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\AvPlusSDevRelMachErr',...
        'w+')
    filenameAvPlusSDevRelMrhoErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\AvPlusSDevRelMrhoErr',...
        'w+')
    filenameAvPlusSDevRelPressErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\AvPlusSDevRelPressErr',...
        'w+')
    
    filenameAvMinusSDevRelTempErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\AvMinusSDevRelTempErr',...
        'w+')
    filenameAvMinusSDevRelMachErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\AvMinusSDevRelMachErr',...
        'w+')
    filenameAvMinusSDevRelMrhoErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\AvMinusSDevRelMrhoErr',...
        'w+')
    filenameAvMinusSDevRelPressErr = open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\AvMinusSDevRelPressErr',...
        'w+')

    # write data to files

    for gridpt = 1:Tot_X_Pts
        fprintf(filenameAvRelTempErr, '#8.6f', AvRelTempErr(gridpt))
        fprintf(filenameAvRelMachErr, '#8.6f', AvRelMachErr(gridpt))
        fprintf(filenameAvRelMrhoErr, '#8.6f', AvRelMrhoErr(gridpt))
        fprintf(filenameAvRelPressErr, '#8.6f', AvRelPressErr(gridpt))
        
        fprintf(filenameAvPlusSDevRelTempErr, '#8.6f', ...
                    AvPlusSDevRelTempErr(gridpt))
        fprintf(filenameAvPlusSDevRelMachErr, '#8.6f', ...
                    AvPlusSDevRelMachErr(gridpt))
        fprintf(filenameAvPlusSDevRelMrhoErr, '#8.6f', ...
                    AvPlusSDevRelMrhoErr(gridpt))
        fprintf(filenameAvPlusSDevRelPressErr, '#8.6f', ...
                    AvPlusSDevRelPressErr(gridpt))
                
        fprintf(filenameAvMinusSDevRelTempErr, '#8.6f', ...
                    AvMinusSDevRelTempErr(gridpt))
        fprintf(filenameAvMinusSDevRelMachErr, '#8.6f', ...
                    AvMinusSDevRelMachErr(gridpt))
        fprintf(filenameAvMinusSDevRelMrhoErr, '#8.6f', ...
                    AvMinusSDevRelMrhoErr(gridpt))
        fprintf(filenameAvMinusSDevRelPressErr, '#8.6f', ...
                    AvMinusSDevRelPressErr(gridpt))
    end

    # close files

    (filenameAvRelTempErr)
    (filenameAvRelMachErr)
    (filenameAvRelMrhoErr)
    (filenameAvRelPressErr)

    (filenameAvPlusSDevRelTempErr)
    (filenameAvPlusSDevRelMachErr)
    (filenameAvPlusSDevRelMrhoErr)
    (filenameAvPlusSDevRelPressErr)

    (filenameAvMinusSDevRelTempErr)
    (filenameAvMinusSDevRelMachErr)
    (filenameAvMinusSDevRelMrhoErr)
    (filenameAvMinusSDevRelPressErr)

    # open files to write residual and its first r time derivatives at throat

    filenameff0 = ...
    open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\ff0Throatvals',...
    'w+')
    filenameff1 = ...
    open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\ff1Throatvals',...
    'w+')
    filenameff2 = ...
    open(...
    'C:\Users\Marek\OneDrive\Documents\ns_q\ff2Throatvals',...
    'w+')

    # write values to files

    for varble = 1:d
        if varble == 1
            for subint = 1:n
                fprintf(filenameff0, '#6.3f',...
                            ff0_throat(varble, subint))
            end
        elseif varble ~= 1
            fprintf(filenameff0, '\n')
            
            for subint = 1:n
                fprintf(filenameff0, '#6.3f', ...
                            ff0_throat(varble, subint))
            end
        end
    end

    for varble = 1:d
        if varble == 1
            for subint = 1:n
                fprintf(filenameff1, '#6.3f',...
                            ff1_throat(varble, subint))
            end
        elseif varble ~= 1
            fprintf(filenameff1, '\n')
            
            for subint = 1:n
                fprintf(filenameff1, '#6.3f', ...
                            ff1_throat(varble, subint))
            end
        end
    end

    for varble = 1:d
        if varble == 1
            for subint = 1:n
                fprintf(filenameff2, '#6.3f',...
                            ff2_throat(varble, subint))
            end
        elseif varble ~= 1
            fprintf(filenameff2, '\n')
            
            for subint = 1:n
                fprintf(filenameff2, '#6.3f', ...
                            ff2_throat(varble, subint))
            end
        end
    end

    # close files and exit

    (filenameff0)
    (filenameff1)
    (filenameff2)
        
    end
    '''
