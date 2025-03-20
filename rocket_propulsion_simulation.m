clc;
clear;

%% Main Script
disp('Enhanced Rocket Propulsion Simulation with Variable O/F Ratio');

% Load configuration and physical constants
config = loadConfiguration();

% Get user inputs
[fuelChoice, oxidizerChoice, ofRatios, thrust, rocketMass, combustionTemp, molarMassExhaust] = getUserInputs(config);

% Validate combustion temperature and molar mass
validateCombustionTemperature(combustionTemp);
validateMolarMass(molarMassExhaust);

% Get fuel and oxidizer properties
fuelProps = getFuelProperties(fuelChoice);
oxidizerProps = getOxidizerProperties(oxidizerChoice);

% Perform calculations
results = calculatePerformance(fuelProps, oxidizerProps, ofRatios, thrust, rocketMass, combustionTemp, molarMassExhaust, config);

% Display results
displayResults(results);

% Plot results
plotResults(results);

%% Functions
function config = loadConfiguration()
    % Load physical constants and configuration
    config.R = 8314; % Universal gas constant in J/(kmol*K)
    config.g0 = 9.81; % Standard gravity in m/s^2
    config.unitConversions.density = 1000; % Convert g/cm^3 to kg/m^3
end

function [fuelChoice, oxidizerChoice, ofRatios, thrust, rocketMass, combustionTemp, molarMassExhaust] = getUserInputs(config)
    % Prompt user for inputs
    disp('Select a fuel:');
    disp('1. Carbon');
    disp('2. Aluminum');
    disp('3. Liquid Methane');
    disp('4. Rocket-Grade Kerosene (RP-1)');
    disp('5. Liquid Hydrogen');
    fuelChoice = input('Enter the number corresponding to your fuel choice: ');
    
    disp('Select an oxidizer:');
    disp('1. Potassium Nitrate (KNO3)');
    disp('2. Ammonium Perchlorate (AP)');
    disp('3. Liquid Oxygen (LOX)');
    oxidizerChoice = input('Enter the number corresponding to your oxidizer choice: ');
    
    ofRatios = input('Enter an array of Oxidizer-to-Fuel (O/F) ratios to analyze (e.g., [2 4 6]): ');
    thrust = input('Enter desired thrust in Newtons (e.g., 1e6 for 1 MN): ');
    rocketMass = input('Enter total rocket mass in kg: ');
    combustionTemp = input('Enter approximate combustion temperature in K (e.g., 2500-3500): ');
    molarMassExhaust = input('Enter approximate molar mass of exhaust in g/mol (e.g., 22-44): ');
end

function validateCombustionTemperature(combustionTemp)
    % Validate combustion temperature range
    if combustionTemp < 2500 || combustionTemp > 4000
        error('Combustion temperature must be between 2500K and 4000K.');
    end
end

function validateMolarMass(molarMassExhaust)
    % Validate molar mass range
    if molarMassExhaust < 15 || molarMassExhaust > 30
        error('Molar mass of exhaust must be between 15 g/mol and 30 g/mol.');
    end
end

function fuelProps = getFuelProperties(fuelChoice)
    % Define fuel properties
    fuels = {'Carbon', 'Aluminum', 'Liquid Methane', 'RP-1', 'Liquid Hydrogen'};
    fuelDensity = [2.267, 2.7, 0.422, 0.82, 0.071]; % g/cm^3
    fuelEnergy = [32.8, 31, 55.5, 43.15, 119.93]; % MJ/kg
    
    % Validate fuel choice
    if fuelChoice < 1 || fuelChoice > length(fuels)
        error('Invalid fuel choice.');
    end
    
    % Return selected fuel properties
    fuelProps.name = fuels{fuelChoice};
    fuelProps.density = fuelDensity(fuelChoice);
    fuelProps.energy = fuelEnergy(fuelChoice);
end

function oxidizerProps = getOxidizerProperties(oxidizerChoice)
    % Define oxidizer properties
    oxidizers = {'KNO3', 'AP', 'LOX'};
    oxidizerDensity = [2.11, 1.95, 1.141]; % g/cm^3
    
    % Validate oxidizer choice
    if oxidizerChoice < 1 || oxidizerChoice > length(oxidizers)
        error('Invalid oxidizer choice.');
    end
    
    % Return selected oxidizer properties
    oxidizerProps.name = oxidizers{oxidizerChoice};
    oxidizerProps.density = oxidizerDensity(oxidizerChoice);
end

function results = calculatePerformance(fuelProps, oxidizerProps, ofRatios, thrust, rocketMass, combustionTemp, molarMassExhaust, config)
    % Initialize results table
    results = table(ofRatios', 'VariableNames', {'OF_Ratio'});
    
    % Calculate mass fractions
    fuelMassFraction = 1 ./ (1 + ofRatios);
    oxidizerMassFraction = ofRatios ./ (1 + ofRatios);
    
    % Calculate mixture properties
    mixtureDensity = (fuelProps.density .* fuelMassFraction + oxidizerProps.density .* oxidizerMassFraction) * config.unitConversions.density; % kg/m^3
    mixtureEnergy = fuelProps.energy .* fuelMassFraction; % MJ/kg
    
    % Estimate specific impulse (Isp) and exhaust velocity
    molarMassExhaust = molarMassExhaust / 1000; % Convert g/mol to kg/mol
    v_e = sqrt((2 * config.R * combustionTemp) / molarMassExhaust); % Exhaust velocity in m/s
    Isp = v_e / config.g0; % Specific impulse in seconds
    
    % Calculate propellant mass flow rate and thrust-to-weight ratio
    propellantMassFlowRate = thrust ./ v_e; % kg/s
    thrustToWeightRatio = thrust ./ (rocketMass * config.g0);
    
    % Store results
    results.MixtureDensity = mixtureDensity';
    results.EnergyContent = mixtureEnergy';
    results.SpecificImpulse = Isp';
    results.ExhaustVelocity = v_e';
    results.PropellantMassFlowRate = propellantMassFlowRate';
    results.ThrustToWeightRatio = thrustToWeightRatio';
end

function displayResults(results)
    % Display results in a formatted table
    disp('Results for varying O/F ratios:');
    disp(results);
end

function plotResults(results)
    % Plot performance parameters vs O/F ratio
    figure;
    subplot(2, 2, 1);
    plot(results.OF_Ratio, results.SpecificImpulse, '-o');
    xlabel('O/F Ratio');
    ylabel('Specific Impulse (s)');
    title('Specific Impulse vs O/F Ratio');
    
    subplot(2, 2, 2);
    plot(results.OF_Ratio, results.ExhaustVelocity, '-o');
    xlabel('O/F Ratio');
    ylabel('Exhaust Velocity (m/s)');
    title('Exhaust Velocity vs O/F Ratio');
    
    subplot(2, 2, 3);
    plot(results.OF_Ratio, results.PropellantMassFlowRate, '-o');
    xlabel('O/F Ratio');
    ylabel('Propellant Mass Flow Rate (kg/s)');
    title('Propellant Mass Flow Rate vs O/F Ratio');
    
    subplot(2, 2, 4);
    plot(results.OF_Ratio, results.ThrustToWeightRatio, '-o');
    xlabel('O/F Ratio');
    ylabel('Thrust-to-Weight Ratio');
    title('Thrust-to-Weight Ratio vs O/F Ratio');
end
