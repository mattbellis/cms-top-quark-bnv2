# Signal: |ΔB| = 2 top decay

Not filled in yet — deliberately. The two background cards are complete and
runnable today; the signal needs the dimension-9 six-fermion UFO
(`t ū d̄ d̄ d̄ d̄`) and the Pythia8 color-flow handling for the ΔB = 2 vertex,
which is the piece that was still open (antijunction topology).

When the model is in place, the card here should follow the same shape as the
backgrounds:

```
import model <BNV_UFO>
define p = g u c d s u~ c~ d~ s~
define j = g u c d s u~ c~ d~ s~
define l- = e- mu-
define vl~ = ve~ vm~

generate     p p > t t~, (t > u~ d~ d~ d~ d~), (t~ > b~ w-, w- > l- vl~)
add process  p p > t t~, (t~ > u d d d d),    (t  > b w+,  w+ > l+ vl)

output signal_bnv -nojpeg
```

Two things worth settling before you generate anything:

1. **The decay chain may not be expressible as a `,`-chain** if the operator
   is a genuine contact interaction with no propagator structure MG5 can
   factorize. If MG5 balks, the fallback is the RAMBO-based LHE decay
   injector: generate `p p > t t~` with one top decayed semileptonically and
   inject the 4-body BNV decay into the LHE afterwards. That also sidesteps
   the `bwcutoff` question for a contact vertex.
2. **The color flow is the real risk, not the kinematics.** Pythia8 has to
   hadronize a color topology with a junction/antijunction. If that is still
   unresolved, you can still get useful *parton-level* kinematics comparisons
   from the LHE files alone while the shower is being sorted — the b-tag
   multiplicity and jet multiplicity arguments do not depend on it.

The two background samples are independent of all of this, so they are the
right thing to validate the pipeline on first.
