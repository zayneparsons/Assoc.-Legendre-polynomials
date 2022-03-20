# Assoc.-Legendre-polynomials
*****The user must have Sympy, Matplotlib, and numpy installed to use this

This is a Python GUI (using Tkinter) that displays the formulas and makes polar plots of the Associated Legendre Polynomials based off the user inputs 
for l and m. Note: this also plots and gives the equations for the regular Legendre Polynomials if the user sets m=0.
  Often times in Physics and Engineering, especially when solving the Laplace equation in spherical coordinates, or solving the Schrodinger equation for 
the Hydrogen atom, the Legendre polynomials and/or the Associated Legendre Polynomials are introduced. Initially, these functions seem very complex 
and foreign. Other transcendental functions such as sine and cosine were also once this way to any Physicist or Engineer, and with time, these functions
become familiar. 
  I intended to make the plots and equations of these formulas easily accessible with this Python script. Instead of looking in a table for these equations
or finding these functions through the Rodrigues formulas, this GUI provides a simple alternative. Additionally, no predefined functions were used in this
script besides some of the more basic functions like cosine. 
  These polynomials are plotted in polar coordinates. The reason for this is that these functions are functions of cosine. I have chosen plots of the 
Associated Legendre Polynomials as functions of cosine of theta because this is the most common form of the polynomials that I have come across in my 
experience as a Physicist. The magnitude of the Associated Legendre Polynomials are taken to be r, and theta ranges from 0 to 2pi. Since the sign of the 
polynomials are not included when we look at the magnitude, the plot displays red whenever the polynomials are positive and blue whenever they are negative.
If things were not done this way (plot the polynomials instead of the magnitude), there would exist a double trace of the same plot. 
  
  
