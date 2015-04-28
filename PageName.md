# Introduction #

The quantref project contains the framework I use to calculate my solutions to the problems I encounter in my coursework at Stevens in their Masters of Financial Engineering program. The intent is to create a viable framework that demonstrates not just how to perform different calculations, but how to structure a multi-lingual valuation framework for derivatives. Error handling is kept to a minimum except where the handling is instructional.


# Details #

Currently, demonstration code is in the file hw2.java. This will be changed as the project formalizes. Valuations are produced by three interacting components:

  1. Valuator - I valuator is a pricing engine that produces values and risk measures given a StateOfTheWorld and an Instrument
  1. StateOfTheWorld - A container for all the market data needed to value an Instrument
  1. Instrument - A representation of a thing for which a price and risk measures can be determined.

Most objects are immutable, or at least _intended_ to be immutable. Perturbations should be done by using a StateOfTheWorld factory to create a PerturbedStateOfTheWorld instance. If you need to add more data to StateOfTheWorld, you should use the factory to create a new StateOfTheWorld that includes your new data.