# High-Pressure Gas Bottle Discharge & Blowdown Simulation

A MATLAB-based numerical simulation and transient thermodynamic modeling toolkit for analyzing blowdown dynamics, mass discharge rates, and pressure-temperature decay during high-pressure gas bottle discharge.

---

## Project Overview
During the rapid blowdown of high-pressure gas bottles, gas expansion induces significant drops in temperature and pressure inside the reservoir. This project models transient thermodynamic behavior, real-gas expansion effects, and nozzle flow dynamics to predict discharge characteristics over time.

### Objectives:
- Model transient pressure (P), temperature (T), and remaining mass (m) over discharge time.
- Simulate choked and subsonic orifice flow regimes.
- Compare numerical simulation output against experimental test data (DFS Test data.xlsx and OFS test data.xlsx).

---

## Repository Structure

- GasBottle.m : Core transient gas discharge solver
- DFS.m : Dynamic flow simulation routine
- OFS200.m : Simulation model for orifice flow test case 200
- OFS260.m : Simulation model for orifice flow test case 260
- DFS Test data.xlsx : Experimental validation data for DFS tests
- OFS test data.xlsx : Experimental validation data for OFS tests
- README.md : Project documentation

---

## Thermodynamic Governing Relations

1. Mass Conservation:
   dm/dt = - m_dot_out

2. Energy Conservation:
   d(m*u)/dt = - (m_dot_out * h_out) + Q_dot_in
   - Accounts for real-gas expansion cooling.
   - Couples ambient convective heat transfer across the vessel boundary.

3. Choked Flow at Orifice:
   m_dot = C_d * A_t * P_0 * sqrt(gamma / (R * T_0)) * (2 / (gamma + 1))^((gamma + 1) / (2 * (gamma - 1)))
   - Transitions to unchoked subsonic expansion once the pressure ratio falls below the critical threshold.

---

## How to Run the Simulations

1. Open MATLAB (R2021a or newer recommended).
2. Set the working folder to the repository directory containing all .m and .xlsx files.
3. Run the primary blowdown solver from the Command Window:
   GasBottle
4. Run specific orifice case studies for test comparisons:
   OFS200
   OFS260
5. Output plots will display the pressure depletion curve, gas cooling profiles, and comparison against measured physical test points.
