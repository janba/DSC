Deformable Simplicial Complex (DSC) method
===

#### The official implementation of the DSC method in 3 dimensions.

The ability to virtually track deformable interfaces has many applications and is therefore of great interest. 
However, Eulerian methods, such as the level set method, tend to suffer from numerical 
difusion and Lagrangian methods are mostly not able to handle topology changes. This was the reason for 
developing a new method, the DSC method. The purpose of the DSC method is to be able to 
model the deformations of objects which are represented by an explicit representation of the interface. 
Furthermore, it is able to model large deformations and handle topology changes.

For **installation** instructions as well as a description of **how to use** the DSC code, see the [DSC wiki](https://github.com/asny/DSC/wiki)

_If you want news about the development and use of the DSC method, please apply for membership of the mailing list at_
https://groups.google.com/forum/#!forum/dsc-development

---
### Online documentation

https://janba.github.io/DSC/

---
### References

For a list of publications, see http://www2.imm.dtu.dk/~janba/DSC-webpage/

For fluid animations, see http://vimeo.com/mkmi/videos

For topology optimization animations, see http://www2.imm.dtu.dk/~asny/Publications.html

---
### Licence

    Deformabel Simplicial Complex (DSC) method
    Copyright (C) 2013  Technical University of Denmark

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    See licence.txt for a copy of the GNU General Public License.
  

The DSC source code includes the external library CGLA which is a part of the [GEL library](http://www2.imm.dtu.dk/projects/GEL/). Permission to distribute the CGLA library as a part of the DSC project has been granted by copyright owner Andreas Bærentzen. For any other purpose the CGLA library is subject to the original license found in CGLA/intro.pdf.

The DSC source code also includes the [SOIL library](http://www.lonesock.net/soil.html) by Jonathan Dummer which is under the MIT license. It is therefore acceptable to include the SOIL package as long as the copyright notice is retained.

Finally, the [TetGen library](http://wias-berlin.de/software/tetgen/) is used for generating the initial simplicial complexes from .obj files. The TetGen library is under the GNU Affero General Public License (a copy of this license is found in the TetGen folder) and is therefore free to distribute.

---
### Contact information

    Asger Nyman Christiansen
    Matematiktorvet, Building 324, room 130
    Technical University of Denmark
    2800 Kgs. Lyngby, Denmark
    e-mail: asny@dtu.dk
