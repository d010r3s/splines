# Splines - 1st AMI's freshman year project 
:skull:
## Table of contents
* [Terms of reference](#terms-of-reference)
* [Cubic spline](#cubic-spline)
* [Parametric cubic spline](#parametric-cubic-spline)
* [Intersection point of two splines](#intersection-point-of-two-splines)
* [Distance between two splines](#distance-between-two-splines)

## Terms of Reference
 
- [x] Implement the construction of cubic spline given the knots and the values at them
- [x] Construct a parametric cubic, allowing one value of x correspond to several values of y
- [x] Implement an algorithm for finding the intersection point of two splines
- [x] Implement an algorithm for finding the minimum distance between two splines

## Cubic Spline
A cubic spline is a spline constructed of piecewise third-order polynomials which pass through a set of m control points. The second derivative of each polynomial is commonly set to zero at the endpoints, since this provides a boundary condition that completes the system of m-2 equations. This produces a so-called "natural" cubic spline and leads to a simple tridiagonal system which can be solved easily to give the coefficients of the polynomials.

![Formula used](https://i.imgur.com/0lU4qyO.png)
![Formula used](https://i.imgur.com/8lzXqOJ.png)
[source](http://statistica.ru/branches-maths/interpolyatsiya-splaynami-teor-osnovy/)


![Cubic Spline](https://blogs.sas.com/content/iml/files/2020/05/cubicInterp1.png)
	
## Parametric Cubic Spline
In the case of parametric cubic
splines, each spline segment is represented by two equations in the independent variable s:
```
x = f1(s) = a_x*(s-s0)^3 + b_x*(s-s0)^2 + c_x*(s-s0) + d_x
x = a_x*t^3 + b_x*t^2 + c_x*t + d_x
t = s - s0

y = f2(s) = a_y*(s-s0)^3 + b_y*(s-s0)^2 + c_y*(s-s0) + d_y
y = a_y*t^3 + b_y*t^2 + c_y*t + d_y
t = s - s0
```
, 

where s0 represents the value of the independent variable s at the beginning of the segment. For convenience,
we’ve made a variable substitution t = s - s0.

[source](https://www.physicsforums.com/attachments/parametric-spline-tutorialv2-pdf.12898/) 


![Parametric cubic spline](https://i.stack.imgur.com/7hbgQ.png)

#### Methods used:
* [Tridiagonal matrix algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)
* [Fitting parametric cubic splines to a set of points](https://www.physicsforums.com/attachments/parametric-spline-tutorialv2-pdf.12898)

## Intersection point of two splines and distance between two splines

#### Methods used:
* Creating an additional third spline merged from the two given
* [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method)
* [A Gradient-Descent Method](https://hal.archives-ouvertes.fr/hal-03854553/file/annpr.pdf)
