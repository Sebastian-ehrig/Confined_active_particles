# Polar active particle on arbitrary curved surfaces

[![Generic badge](https://img.shields.io/badge/arXiv-1610.05987-green.svg)](https://arxiv.org/abs/1610.05987)
[![Generic badge](https://img.shields.io/badge/Phys.Rev.E-95.062609-yellow.svg)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.062609)
![Generic badge](https://img.shields.io/badge/Matlab-R2015b-blue.svg)
[![Generic badge](https://img.shields.io/badge/License-AGPL3.0-orange.svg)](https://github.com/Sebastian-ehrig/Confined_active_particles/blob/main/LICENSE)

---

The aim of the model is to explores the effect of surface curvature on the collective behaviour of active particles.

It is part of a wider research initiative in the context of pattern formation in biological systems.

The results of the simulation for ellipsoidal surfaces have been published in ["Curvature-controlled defect dynamics in active systems"]( https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.062609).

The accepted version of the paper is available for free on arxiv: https://arxiv.org/abs/1610.05987.

---

## Theoretical Model

The model described in the publication is a Vicsek type model (Vicsek et al. 1995, Physical review letters 75(6): 1226) of spherical active particles with a fixed radius confined to the surface of an ellipsoid. Particle interactions are modelled through forces between neighbouring particles that tend to align their velocities (adapted from Szabo et al. 2006, Physical Review E 74(6): 061908).

The particle motion on the curved surface is performed by an unconstrained motion in the tangential plane followed by a projection onto the surface. To be able to use the model on arbitrary surfaces, surfaces are approximated with triangulated meshes.

All simulations have been performed in Matlab R2015b by solving an overdamped differential equations of motion using a forward Euler integration method with a fixed time step. Model paprameters have been chosen such that the study is in the regime of low noise and low energy and particles interact virually as hard spheres. 

For more details of the model please see publication.

---
### Example Videos

Supplemental material Movies of the paper ["Curvature-controlled defect dynamics in active systems"](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.062609).
Points of constant normal curvature (umbilics) are highlighted as red spheres, black arrows are the particle orientation vectors and the surface is color-coded according to the vortex order parameter (VOP, dark blue for VOP=0 and yellow for VOP=1).

##### Movie 1

Example of particle and vortex dynamics on a prolate spheroid with aspect ratio x/z = 0.25.

<a href="https://www.youtube.com/watch?v=M6JfDUZr97I" target="_blank">
 <img src="http://img.youtube.com/vi/M6JfDUZr97I/maxresdefault.jpg" alt="Watch the video" width="480" height="270" border="0" />
</a>

<center>
 <iframe width="560" 
         height="315" 
         src="https://www.youtube.com/watch?v=M6JfDUZr97I" 
         title="YouTube video player" 
         frameborder="0" 
         allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" 
         allowfullscreen>
 </iframe>
</center>

---

## Reference

Ehrig, S., Ferracci, J., Weinkamer, R., & Dunlop, J. W. Curvature-controlled defect dynamics in active systems. Phys. Rev. E 95, 062609 (2017)
https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.062609

Ehrig, S., Ferracci, J., Weinkamer, R., & Dunlop, J. W. (2016). Curvature-controlled defect dynamics in active systems. arXiv preprint arXiv:1610.05987.
https://arxiv.org/abs/1610.05987
