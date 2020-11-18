# qsimulatR

A simple quantum computer simulator in R.

## Installation

Install the [programming language R](https://www.r-project.org/) if you have not done so yet.

Download the source code, e.g. with
```
git clone https://github.com/HISKP-LQCD/qsimulatR.git
```
and go into the newly created directory `qsimulatR`.
Execute:
```
./install
```
You might have to install additional packages. Just use `install.packages()` for any packages recommended in a possible error message.

The `install` script might not work under Windows. To circumvent this, you can download the source package from `github`. Then
```
install.packages("qsimulatR.zip", repos=NULL, type="source")
```

The library `qsimulatR` is now available in your R installation and can be loaded with:
```
library(qsimulatR)
```

Check for updates regularly. Do so by going into the directory `qsimulatR` and executing
```
git pull
./install
```

## Usage
A detailed usage description can be found in `qsimulatR.pdf`. We provide many useful examples in the R Markdown format in `vignettes`.

Your first very simple program with `qsimulatR` might look like this:
```
library(qsimulatR)

# generate a quantum state with 2 qubits, initialised to |00>
x = qstate(nbits=2)
# display the state
x

# apply the Hadamard gate to the first (right) qubit
y = H(1) * x
y

# apply a controlled NOT
z = CNOT(c(1,2)) * y
z

# draw the resulting circuit
plot(z)

# project onto a single compute basis state
res = measure(z)
# draw the circuit
plot(res$psi)

# perform the measurement many times and plot the outcome
dist = measure(z, rep=1000)
hist(dist)
```
