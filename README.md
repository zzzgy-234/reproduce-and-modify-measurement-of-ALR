# CEPC Polarimetry and Weak Mixing Angle Measurement

This repository contains research and simulations focused on determining the effective weak mixing angle, \(\sin^2\theta_{eff}\), with unprecedented precision at the GigaZ option of the CEPC (Circular Electron Positron Collider). Achieving this precision requires keeping the beam polarization uncertainty below 0.1%. This level of precision can be achieved through a modified [Blondel scheme](https://www.sciencedirect.com/science/article/abs/pii/0370269388908696), as Compton polarimetry alone is insufficient.

## Key Highlights:

- **Blondel Scheme**: The Blondel scheme uses electron-positron helicity combinations to measure the scattering cross section, avoiding the need for direct polarization measurement. However, at LEP, a single-ring collider, both electrons and positrons were polarized simultaneously, limiting the available helicity combinations. Additionally, Blondel’s original assumptions about polarization were idealized.

- **CEPC Double-Ring Collider**: We reproduced the Blondel scheme and leveraged the CEPC’s double-ring collider design to provide a richer set of electron-positron helicity combinations, allowing for more precise measurements.

- **Polarization Correction**: An additional correction condition for polarization was introduced to improve the consistency of the original Blondel scheme with real-world conditions in polarization asymmetry measurements.

- **Monte Carlo Simulations**: Based on Monte Carlo algorithms, the scattering cross sections of electron-positron collisions were computed with varying luminosities and helicities.

- **Genetic Algorithm Optimization**: A genetic algorithm was used for real-time optimization of the luminosity distribution, ensuring measurement precision.

- **Error Estimation**: With electron polarization \( P^- = 70\% \) and positron polarization \( P^+ = 30\% \), along with the provision of over ten thousand bunches at CEPC, the absolute error of \( A_{LR} \) was limited to \( 10^{-5} \). This assumes a 5% relative error in polarization measurement and no deviation from the expected polarization of 0.05.

---

This repository contains the code, data, and simulations used in this study. Contributions and discussions are welcome.
