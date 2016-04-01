Playing around with scientific calculation in different programming languages.

These tests are not really benchmarks, but can give a rough feeling about what orders we are talking.

Fortran: 150ms
Java   : 460ms
Python : 12500ms

One of the surprising parts was, how unreadable the Java-Code was in comparison to Python (I ported from python
to Java). Python on the other hand had the problem that I really had to step through the code with a debugger
to figure out which datatypes the functions used.

More material on Java in HPC:

High-performance computing in Java: the data 
processing of Gaia
X. Luri & J. Torra ICCUB/IEEC
http://www.spscicomp.org/ScicomP15/slides/astro/torra.pdf
Paper: http://www.aspbooks.org/publications/434/135.pdf


"Current State of Java for HPC"
Brian Amedro 1 Vladimir Bodnartchouk 2 Denis Caromel 1 Christian Delbe 1 Fabrice Huet 1 Guillermo Taboada 3, * 
https://hal.inria.fr/inria-00312039/document

http://www.des.udc.es/~gltaboada/papers/JAVAHPC-CIEMAT-2011.pdf



Todo
====

Do it in scala
http://www.scalanlp.org/