# AE-Senior-Design
Senior Design Code for Project Goose
This is a repository for the Saint Louis University 2022-2023 Senior Design project: Project Goose
The codes present are the basis for our aerodynamic calculations

## Setting up Git in MATLAB
1. In an **EMPTY** MATLAB folder, right click inside the folder browser 
2. Navigate to "source control" and click "manage files" (seen in picture below)
![Step1](https://imgur.com/skGyjZP)
3. On the Github page you are reading this on, click the green code button and copy the link that is displayed (the HTTPS one)
![Step2](https://imgur.com/9lW3jcB)
4. Paste the link into "repository path" in the MATLAB window seen below
![Step3](https://imgur.com/FeGlWyR)
5. Click "Retrieve

## Drag Code
Drag.m is the code for calculation of aircraft drag coefficient, lift coefficient, takeoff and landing distances, and cruise lift and drag.
It is dependant on

* DragBuildup.m
- PlaneInfo.mat

To run DragBuildup.m in a code, input 
```
[CL, CD] = DragBuildup("Angle of Attack in degree", PlaneInfo.mat)
```

## Constraint Analysis 
ElectricConstraint.m runs a constraint analysis for an electric aircraft. It also calculates a very rough weight analysis. 

## Stabilizer Sizing
RoughSizing.m calculates the size of the horizontal and vertical stabilizers based off volume coefficients.

