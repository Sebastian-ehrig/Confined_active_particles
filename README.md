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

Markdown style

[![Watch the video](https://img.youtube.com/vi/nTQUwghvy5Q/default.jpg)](https://youtu.be/nTQUwghvy5Q)

HTML style

<a href="http://www.youtube.com/watch?feature=player_embedded&v=nTQUwghvy5Q" target="_blank">
 <img src="http://img.youtube.com/vi/nTQUwghvy5Q/mqdefault.jpg" alt="Watch the video" width="240" height="180" border="10" />
</a>

Thumbnail

![alt text][image]
[image]: https://img.youtube.com/vi/M6JfDUZr97I/1.jpg "alt_txt"

https://www.youtube.com/watch?v=M6JfDUZr97I

---

## Reference

Ehrig, S., Ferracci, J., Weinkamer, R., & Dunlop, J. W. Curvature-controlled defect dynamics in active systems. Phys. Rev. E 95, 062609 (2017)
https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.062609

Ehrig, S., Ferracci, J., Weinkamer, R., & Dunlop, J. W. (2016). Curvature-controlled defect dynamics in active systems. arXiv preprint arXiv:1610.05987.
https://arxiv.org/abs/1610.05987
