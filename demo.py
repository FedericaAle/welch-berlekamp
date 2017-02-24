from copy import deepcopy
from welchberlekamp import makeEncoderDecoder
from finitefield.finitefield import FiniteField

## Setup ##
n = 3   # Number of Packets we want to send
k = 1   # Maximum number of errors we can tolerate
p = 11   # We're working on finite field modulo p
Fp = FiniteField(p)
enc, dec, solveSystem = makeEncoderDecoder(n + 2 * k, n, p)

## Encode Message like this ##
message = [3, 4, 0]
encoded = enc(message)
print("The encoded message is:")
print(encoded)

## Decode Message like this ##
received = deepcopy(encoded)
received[0]=Fp(3) #Changing a message in received.
received[1]=Fp(7) #Changing a message in received.
received[2]=Fp(0) #Changing a message in received.
received[3]=Fp(2) #Changing a message in received.
received[4]=Fp(10) #Changing a message in received.
print("The received message is:")
print(received)
dec(received)