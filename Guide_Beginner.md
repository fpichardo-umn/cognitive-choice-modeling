# Computational Modeling Guide for Newcomers

*Understanding the fundamentals of computational cognitive modeling*

---

## What is Computational Modeling?

**The Core Idea:** Instead of just describing what people do, we build mathematical models that explain *how* they do it.

### A Real Example: The Iowa Gambling Task

Imagine you're studying how people learn to make better decisions. In the Iowa Gambling Task (IGT), participants repeatedly choose from four decks of cards:

- **Deck A:** Frequent large wins (+100), frequent large losses (-250; 80%), net loss
- **Deck B:** Frequent lage wins (+100), rare huge losses (-250; 10%), net loss
- **Deck C:** Frequent small wins (+50), frequent small losses (-50; 80%), net gain
- **Deck D:** Frequent small wins (+50), rare small losses (-250; 10%), net gain

**The Pattern:** Most people start by exploring all decks, but gradually learn to prefer Decks C and D (the "good" decks) over Decks A and B (the "bad" decks).

**The Question:** What mental processes drive this learning?

---

## Why Traditional Analysis Isn't Enough

### Traditional Approach:
```
Average performance in last 20 trials: 75% good deck choices
Conclusion: "People learned to prefer good decks"
```

**Problems with this approach:**
- **No insight into process:** We know *what* happened, not *how*
- **Individual differences ignored:** Some people learn fast, others slow, others never learn
- **No theory testing:** We can't distinguish between competing explanations
- **No prediction:** We can't predict behavior in new situations

---

## The Computational Modeling Solution

Instead of just measuring outcomes, we propose **specific theories about the mental processes** and test them quantitatively.

### Example Theories (Models):

**Theory 1: Expected Value Model**
> "People track the average outcome of each deck and choose accordingly"

**Theory 2: Prospect Valuation Learning Model**  
> "People weight losses more heavily than gains, and have subjective value functions"

**Theory 3: Outcome Representation Learning Model**
> "People have separate learning systems for rewards and punishments"

### How It Works:

1. **Translate theories into math:** Each theory becomes a computational model with specific equations
2. **Fit models to data:** Find parameter values that best reproduce each person's choices
3. **Compare models:** Which theory best explains the data across participants?
4. **Interpret parameters:** What do the best-fitting values tell us about individuals?

---

## Understanding Model Components

### The Data: What We Observe
```
subjID trial  choice  wins  losses
100    1      2       100    0
100    2      1       100    0  
100    3      1       250    50
100    4      4       50     0
100    5      3       50     50
...    ...          ...      ...
```

### The Model: What We Don't Observe

**Hidden mental processes that might explain the choices:**

```
Trial 1: Person sees all decks as equal → explores randomly
Trial 3: Big loss from Deck 1 → updates expectation downward  
Trial 10: Deck 1 expectation now negative → avoids it
Trial 50: Clear preference for Decks 3 & 4 → exploits knowledge
```

**The model's job:** Predict choice probabilities based on these hidden processes.

---

## Available Model Families

### Reinforcement Learning (RL) Models

**Core Idea:** People learn from experience by updating their expectations.

#### **Expected Value (EV) Models**
- **What it assumes:** People track the average outcome of each option
- **Learning process:** After each choice, update the expected value toward the observed outcome
- **Key parameters:**
  - `update`: How quickly expectations change (0 = no learning, 1 = complete updating)
  - `con`: How consistently people choose the best option (low = random, high = deterministic)
  - `wgt_rew`: Weight given to gains
  - `wgt_pun`: Weight given to losses

#### **Prospect Valuation Learning (PVL) Models**
- **What it adds:** Psychological value functions from prospect theory
- **Key insight:** A $100 gain doesn't feel the same as a $100 loss
- **Additional parameters:**
  - `A`: Choice consistency (similar to `con` but different scale)
  - `w`: Loss aversion (how much worse losses feel than equivalent gains)
  - `a`: Risk attitude (concave vs. convex value function)

#### **Value-Plus-Perseveration (VPP) Models** 
- **What it adds:** Tendency to repeat recent choices regardless of outcomes
- **Key insight:** People sometimes stick with options due to habit/inertia
- **Additional parameters:**
  - `ep`: Perseveration strength
  - `pattern_weight`: Contribution of perseveration vs. value

#### **Outcome Representation Learning (ORL) Models**
- **What it assumes:** Separate learning systems for rewards and punishments
- **Key insight:** Brain systems for "seeking rewards" vs. "avoiding losses" work differently
- **Different parameters:**
  - `Arew`: Learning rate for rewards
  - `Apun`: Learning rate for punishments  
  - `betaF`: Choice consistency for frequency-based learning
  - `betaP`: Choice consistency for payoff-based learning

#### **Value and Sequential Exploration (VSE) Models**
- **What it adds:** Explicit exploration strategies beyond random choice
- **Key insight:** People sometimes deliberately explore to gather information

#### **Learning Rules Within Models:**
- **Delta learning:** Update expectations toward most recent outcome
- **Decay learning:** Expectations gradually fade toward baseline over time
- **Both:** Combination of delta learning + decay

### Sequential Sampling Models (SSM)

**Core Idea:** Decision-making is a noisy process that unfolds over time, explaining both choices and reaction times.

#### **Drift Diffusion Model (DDM)**
- **What it assumes:** Evidence accumulates toward decision boundaries
- **Key insight:** Both choice and reaction time emerge from the same underlying process
- **Parameters:**
  - `boundary`: How much evidence needed before deciding (higher = slower, more accurate)
  - `drift`: Rate of evidence accumulation toward correct choice
  - `tau`: Non-decision time (perception, motor response)

### Hybrid Models (RL + SSM)

**Core Idea:** Combine learning from experience with realistic decision timing.

- **Example:** `pvldelta_ddm` combines PVL learning with DDM decision process
- **What you get:** Predictions for both choice patterns AND reaction time distributions
- **Current status:** Well-developed for modified IGT, under development for standard IGT

---

## Model Parameters: What They Mean

### Learning Parameters
```
update = 0.1    # Slow learner: gradual belief updating
update = 0.9    # Fast learner: rapid belief updating
```

### Motivational Parameters  
```
wgt_rew = 0.8   # Moderately motivated by gains
wgt_pun = 0.9   # Strongly motivated to avoid losses (loss averse)

# Or in PVL models:
w = 2.0         # Losses feel twice as bad as equivalent gains
a = 0.8         # Risk averse (concave value function)
```

### Control Parameters
```
con = 3.0       # Very consistent: almost always chooses best option
con = 0.5       # Somewhat random: often chooses suboptimal options
con = -1.0      # Very random: choice barely relates to values
```

### Decision Time Parameters (DDM models)
```
boundary = 1.5  # Cautious: collects lots of evidence before deciding
boundary = 0.8  # Impulsive: decides quickly with less evidence
drift = 2.0     # Good at discriminating: strong evidence accumulation
tau = 0.3       # 300ms non-decision time (perception + motor response)
```

---

## The Research Process

### 1. **Model Selection**
- Choose models that represent theories you want to test
- Consider what psychological processes you're interested in
- Start simple (EV model) before moving to complex (hybrid models)

### 2. **Model Fitting** 
- Use statistical methods (Bayesian estimation) to find parameter values
- Each person gets their own parameter estimates
- Good fits require proper convergence and adequate sampling

### 3. **Model Validation**
- **Parameter Recovery:** Can the model reliably estimate its own parameters?
- **Posterior Predictive Checks:** Does the fitted model produce realistic behavior?
- Both steps are crucial before trusting results

### 4. **Model Comparison**
- Compare multiple models using information criteria (LOOIC/WAIC)
- Consider both statistical fit and psychological interpretability
- Best statistical fit doesn't always mean best theory

### 5. **Interpretation**
- What do the parameter values tell us about individuals?
- How do parameters relate to other measures (personality, brain activity, etc.)?
- Can we predict behavior in new situations?

---

## Key Insights from Computational Modeling

### Individual Differences
Instead of: *"People learned to prefer good decks"*
We get: *"Fast learners (high update) showed early preference shifts, while slow learners (low update) required more experience. Loss-averse individuals (high wgt_pun) showed stronger avoidance of bad decks."*

### Process Understanding  
Instead of: *"Performance improved over time"*
We get: *"Choice consistency (con parameter) increased across trials, while learning rates (update) remained stable, suggesting improved decision-making rather than faster learning."*

### Theory Testing
Instead of: *"Learning occurred"*  
We get: *"The PVL model outperformed the EV model (ΔLOOIC = 15.3), supporting prospect theory predictions about asymmetric gain/loss processing."*

### Clinical Applications
Instead of: *"Patients performed worse"*
We get: *"Patients showed intact learning rates but reduced choice consistency, suggesting decision implementation deficits rather than learning impairments."*

---

## Getting Started: Key Concepts

### Essential Understanding:
1. **Models are theories:** Each represents a specific hypothesis about mental processes
2. **Parameters are individual differences:** They capture how processes vary across people  
3. **Validation is crucial:** Never trust a model without testing it thoroughly
4. **Comparison is key:** Multiple models provide stronger inference than single models
5. **Interpretation matters:** Statistical fit alone doesn't guarantee psychological insight

### Common Pitfalls to Avoid:
- **Skipping validation:** Always test parameter recovery and behavioral realism
- **Over-interpreting parameters:** Ensure the model can actually estimate them reliably
- **Ignoring model comparison:** Single models can be misleading
- **Forgetting psychology:** Statistical fit doesn't guarantee theoretical insight
- **Rushing to complexity:** Master simple models before attempting hybrid approaches

### Next Steps:
Once you understand these concepts, you'll be ready to learn the practical aspects of fitting models to data, validating results, and interpreting findings. The computational modeling pipeline provides the tools to implement these concepts rigorously and efficiently.

---

*This guide provides the conceptual foundation for computational cognitive modeling. The next step is learning to use the modeling pipeline to implement these concepts with real data.*
