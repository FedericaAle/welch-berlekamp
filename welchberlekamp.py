# an encoder and decoder for Reed-Solomon codes with coefficients in Z/p for a prime p
# decoder uses the Berlekamp-Welch algorithm

# for solving a linear system
from linearsolver.linearsolver import someSolution

from finitefield.finitefield import FiniteField
from finitefield.polynomial import polynomialsOver


def makeEncoderDecoder(n, k, p):
   if not k <= n <= p:
      raise Exception("Must have k <= n <= p but instead had (n,k,p) == (%r, %r, %r)" % (n,k,p))

   Fp = FiniteField(p)
   Poly = polynomialsOver(Fp)
   maxE = ((n - k) // 2)  # maximum allowed number of errors

   # message is a list of integers at most p
   def encode(message):
      if not all(x < p for x in message):
         raise Exception("Message is improperly encoded as integers < p. It was:\n%r" % message)

      def row(i, b):
         return [Fp(i ** (j)) for j in range(k)] + [Fp(b)]
      system = [row(i, message[i]) for i in range(k)]
    
      interpolated = someSolution(system, freeVariableValue=1)
      thePoly = Poly(interpolated)
      print("The polynomial encoding the message is:")
      print(thePoly)
      return [thePoly(Fp(i)) for i in range(n)]


   def solveSystem(encodedMessage, debug=True):
      for e in range(maxE, 0, -1):
         ENumVars = e+1
         QNumVars = e+k
         def row(i, a, b):
            return ([b * a**j for j in range(ENumVars)] +
                    [-1 * a**j for j in range(QNumVars)] +
                    [0]) # the "extended" part of the linear system

         system = ([row(i, a, b) for (i, (a,b)) in enumerate(encodedMessage)] +
                   [[Fp(0)] * (ENumVars-1) + [Fp(1)] + [Fp(0)] * (QNumVars) + [Fp(1)]])
                     # ensure coefficient of x^e in E(x) is 1

         if debug:
            print("\nsystem is:\n\n")
            for row in system:
               print("\t%r" % (row,))

         solution = someSolution(system, freeVariableValue=1)
         E = Poly([solution[j] for j in range(e + 1)])
         Q = Poly([solution[j] for j in range(e + 1, len(solution))])

         if debug:
            print("\nreduced system is:\n\n")
            for row in system:
               print("\t%r" % (row,))

            print("solution is %r" % (solution,))
            print("Q is %r" % (Q,))
            print("E is %r" % (E,))

         P, remainder = Q.__divmod__(E)
         if remainder == 0:
            return Q, E

      raise Exception("found no divisors!")


   def decode(encodedMessage):
      encodedMessage = [[Fp(i), encodedMessage[i]] for i in range(len(encodedMessage))]
      Q,E = solveSystem(encodedMessage)

      Pcoefs, remainder = Q.__divmod__(E)
      if remainder != 0:
         raise Exception("Q is not divisibly by E!")
      P = Poly(Pcoefs)
      return [P(Fp(i)) for i in range(k)]
      


   return encode, decode, solveSystem
