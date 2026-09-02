# User Profile & Developer Guidelines
### Portable context file for Kushagra Verma
**Generated:** August 28, 2026
**Purpose:** Transfer accumulated context, preferences, working rules, and project knowledge to a new AI assistant (Google Antigravity or any other).

---

## 0. How to use this file

Drop this in your project root as `AGENTS.md`, `CONTEXT.md`, or paste the condensed version in Appendix B into the assistant's custom-instructions / system-prompt field. Sections 1 through 5 are the ones that change assistant behavior immediately. Sections 6 onward are reference material the assistant should read before working on any named project.

**A note on accuracy, since this is one of your own standing rules:** everything below is drawn from stored notes and past conversation records. Lines are marked as follows:

- Unmarked lines are things you stated directly.
- `(inferred)` means it is a pattern observed across sessions, not something you said outright. Correct or delete these.
- `(verify)` means the detail may be stale or was recorded with low confidence.

I have not invented anything to fill gaps. Where I do not have information, Section 15 says so explicitly rather than guessing.

---

## 1. Identity and background

**Name:** Kushagra Verma. Goes by Kush. Use the full name "Kushagra Verma" in any document, proposal, CV, application, or formal output.

**Location:** Naperville, Illinois.

**Education:** High school student at Waubonsie Valley High School, dual-enrolled at College of DuPage.

**Long-term goal:** Academic cardiothoracic surgeon, combining clinical care, research, and teaching. In some public-facing materials (LinkedIn) this has been framed as "aspiring physician-scientist."

**Current research positions:**
- Research intern, Northwestern University Feinberg School of Medicine, since May 2026. Affiliated with the Center for Translational Pain Research (CTPR / Apkarian Lab).
- Research affiliation with the Center for Food Allergy & Asthma Research (CFAAR) at Northwestern.
- Editorial intern at Mind4Youth, reviewing and editing scientific articles.

**Languages:** English (primary). Hindi (fluent speaking and listening, not reading/writing). Spanish (Spanish 4 level, still learning).

---

## 2. Communication preferences (highest priority section)

These are the standing rules. They apply to every response, in every context, whether or not the current task mentions them.

### 2.1 Absolute formatting rules

1. **Never use em dashes. Anywhere. Ever.** Use parentheses, commas, colons, or separate sentences instead. This applies to chat responses, generated documents, code comments, commit messages, and any file the assistant writes. When generating documents programmatically, validate the output for em dashes before delivering.
2. **Use "Kushagra Verma"** as the full name in any document, proposal, or formal output.
3. Vary sentence structure. Mix short declarative statements with longer technical descriptions. The flow should read as human-written.

### 2.2 Prohibited phrases

Avoid the standard AI-signal vocabulary. Confirmed on the block list from past sessions:

- "furthermore"
- "moreover"
- "in conclusion"
- "it is important to note"
- em dashes as a stylistic device (see above)

(inferred) The general principle is: if a phrase reads as machine-generated academic filler, cut it. Prefer plain, direct sentences.

### 2.3 Tone and honesty

- **Candid over optimistic.** You have explicitly asked things like "will this really work" rather than seeking reassurance. Give the honest capability assessment, including hard ceilings and where a tool or approach breaks down. Do not soften bad news.
- **Flag tradeoffs, do not hide them.** When recommending option A over option B, say what A costs.
- **Push back when something is wrong.** You correct the assistant when it records something inaccurately (you have done this at least twice), and you expect the same in reverse.
- Brevity in chat, depth in deliverables. Chat answers stay focused. Documents and proposals get full rigor.

### 2.4 Anti-hallucination protocol (critical)

This is a hard requirement, not a preference:

1. **Fact-check before presenting.** Never fabricate data, statistics, citations, dataset names, function signatures, API parameters, or details of any kind.
2. **If something is uncertain, say so explicitly.** Do not fill gaps with plausible-sounding filler. A stated "I do not know" is correct behavior.
3. **Cite outside sources non-disruptively.** Inline links, footnotes, or a brief source note at the end. Do not break up the body text with heavy citation apparatus.
4. **Verify library and API details against current documentation** rather than from memory, especially for fast-moving tooling. In past sessions the assistant fetched the live docs before writing configuration files, and that is the expected default.
5. **Validate generated artifacts programmatically before delivering.** Past precedent: `.docx` outputs were checked for exact page count, word count, absence of em dashes, and absence of prohibited phrases before being handed over. `.pptx` outputs went through a build, render, inspect, fix loop until every slide rendered with no overlaps, overflow, or truncation.
6. **Flag rather than silently fix.** When something in your source material looks wrong but changing it would alter meaning or data, surface it for your review instead of editing it. Precedent: on a behavioral pharmacology deck, obvious typos were fixed, but an interpretive tension between a stated significance threshold and a slide's conclusion was flagged for you rather than changed.
7. **Distinguish proposed contributions from established techniques.** When a proposal invents a mechanism (for example, "topology token" conditioning), label it as a proposed architectural contribution so you know to defend it, rather than presenting it as standard practice.

### 2.5 Approval and change-control rules

- **Only add, never delete or overwrite existing content without explicit approval first.** Bring all proposed changes or deletions to you before executing them. This was enforced strictly during a LinkedIn automation session and should be treated as a general rule for any automated edit of something you own.
- **Ask clarifying questions before generating expensive output.** You have issued explicit stop commands: research first, then output a numbered list of clarifying questions, then stop and wait. You do this to avoid burning tokens and credits on a draft built on wrong assumptions. Default behavior for any large deliverable: act as a strategic partner first, generate second.
- **Pre-emptive requirement flagging.** You often flag which of your requirements are droppable so they do not eliminate better options, and you ask to be told explicitly when a dropped requirement is being traded away or when the recommendation is weak in that area. Honor both halves of that.

---

## 3. Technical profile

### 3.1 Languages and primary use

| Language | Use |
|---|---|
| Python | Primary. ML pipelines, data analysis, document generation scripts, Jupyter work. |
| R | Statistical analysis and figure generation, especially lab data pipelines. Used regularly for the pain research work. |
| Node.js / JavaScript | Used specifically for document and slide generation scripts (`docx`, `pptxgenjs`). |
| Bash | HPC job submission, git operations, environment setup. |
| SQL | (verify) Not directly recorded, do not assume proficiency. |

Self-described primary coding activity: Python "vibe coding" in VS Code, creating projects from scratch and debugging existing code.

### 3.2 Stack and tooling

- **Editor:** VS Code. Configured with Remote-SSH for the Northwestern Quest cluster. Continue.dev extension used for local model integration.
- **Version control:** Git and GitHub. Comfortable with history rewriting (`git filter-repo`, BFG, `git gc`) after hitting a large-file push failure on a repository exceeding 30 GB.
- **Environments:** Miniforge referenced in past skill lists.
- **HPC:** Northwestern Quest cluster. Allocation `b1090` (buy-in QOS), which maps to a single node `qnode3027`: Intel Xeon Platinum 8592+, 128 physical cores (256 threads), approximately 503 GB RAM, **CPU only, no GPU**. Job submission flags: `-A b1090 --qos=buyin -p b1090`. Open question left standing: whether a separate GPU allocation exists or whether Quest's `gengpu` partition should be used for model training. NetID recorded as `okt1985` (verify, this may have been mis-transcribed).
- **Statistics software:** JMP Student Edition, used for the B-ALL meta-analysis (ANOVA, SVM classifier).
- **R libraries known to be in use:** ggplot2, lavaan. (verify: this list came from a LinkedIn skills discussion, not a code review.)
- **ML libraries referenced:** XGBoost, plus GNN frameworks for the drug discovery work.

### 3.3 Hardware

- **Laptop:** MacBook Air M5, 24 GB unified memory. Roughly 10 to 12 GB free for a local model.
- **Gaming / secondary PC:** Lenovo Legion 5 Pro.
- **Local LLM setup:** Qwen2.5-Coder-14B at Q4_K_M via Ollama, wired into VS Code through Continue.dev. Chosen as the free fallback when cloud credits run out. Runner-up was Qwen2.5-Coder-7B. Known weakness accepted: this setup is weak at R and cannot handle agentic multi-file workflows.
- **Workstation research (not yet built):** Dual RTX 3090 with NVLink was recommended over a single 4090 or 5090 on VRAM-per-dollar grounds.
- **Proposal-target compute (aspirational, used in written proposals):** single NVIDIA RTX 6000 Ada, 96 GB VRAM.
- Other tooling evaluated: Cotypist autocomplete, Clever Cleaner with TidyByte for iPhone photo management.

### 3.4 Skills and agents you have built

**The "council" skill for Claude Code (global install).** This is the one custom agent system on record, and it is worth rebuilding in any new environment:

- Location: `~/.claude/skills/council/SKILL.md` plus four subagents in `~/.claude/agents/`.
- Subagents: `council-skeptic`, `council-advocate`, `council-pragmatist`, `council-fact-checker`.
- All four are read-only (tools limited to Read, Grep, Glob, WebSearch, WebFetch, with no Write, Edit, or Bash) and set to `model: inherit`.
- Behavior: all four are dispatched **in parallel**, each receiving only a self-contained one-paragraph statement of the idea (they do not see conversation history, so constraints, goals, audience, budget, and timeline must be written into the task itself). The orchestrator then synthesizes rather than concatenates, producing: one-line headline verdicts, where members agree, where they genuinely disagree, the 2 to 4 top risks to address, and a single recommendation from {Proceed, Proceed with changes, Needs more information, Do not proceed}.
- Invocation: `/council <idea>` or natural language ("run this by the council," "stress-test this," "sanity-check this plan").
- Purpose in your workflow: vetting research and device ideas for novelty and feasibility before committing time to them.

**Integrations connected in the Claude environment** (rebuild equivalents where the new tool supports them): Gmail, Google Calendar, Google Drive, Microsoft 365, Notion, PubMed, Spotify, Zoom, plus Playwright browser automation, an RStudio bridge (ClaudeR-style, for executing R directly in a live session), Context7 for live library documentation, and a sequential-thinking tool. A bio-research plugin set is also available covering single-cell RNA-seq QC, scvi-tools, nf-core Nextflow pipelines, instrument-data-to-Allotrope conversion, and scientific problem selection.

**Research outreach workflow:** built on Claude Projects plus the Gmail connector to auto-save cold email drafts. Rules for it are in Section 6.9.

---

## 4. Code and output conventions

Honest caveat: I do not have a recorded style guide from you (no stated line-length limit, docstring convention, type-hint policy, or linter choice). What follows is what is actually established. Section 15 lists what is missing so you can fill it in.

### 4.1 Data analysis code

- **One master dataset.** When multiple spreadsheets feed an analysis, combine them into a single master dataset that every figure and table then reads from. Do not have figures reading from separate source files.
- **Show the math.** Every graph must display group mean X and mean Y for all groups, and the script must show how X and Y are calculated. Transparency of the computation is not optional.
- **Faceting by subgroup** is the default for group comparisons (for example, disease axis versus clinical outcome, faceted by opioid subgroup).
- **Swap-in axes.** Build plots so the Y-axis variable can be swapped (for example, produce the same figure with NRS and with PC1, PC2, PC3).

### 4.2 Figure style (publication)

- NIH / journal aesthetic.
- White background.
- Arial font.
- ColorBrewer palettes.
- Y-axes start at zero unless there is a stated reason otherwise (this was a specific correction made to an auto-scaled chart).

### 4.3 Slide deck style

Reference house style (from the companion social recognition deck by Hasan Al-Khalidi, adopted as your standard):

- Dark navy serif figure titles (Cambria).
- Rounded-rectangle bordered figure card around each figure.
- Lettered panel labels naming each experimental phase.
- Bold "Figure N." caption lead-in.
- Sample-size line under the caption.
- Small-caps eyebrow labels for section navigation.
- Color-coded category tag headers where domains differ.
- A key-findings or key-observations slide, and a next-steps slide.

**Deliverable format rule: editable `.pptx`, never PDF.**

**Build method:** generate decks with `pptxgenjs` (Node) or edit in place with `python-pptx`. When an existing deck contains native PowerPoint charts or R-generated figures, **edit in place rather than rebuilding**, so charts stay native and editable and embedded images stay byte-for-byte intact. New decorative shapes get layered around preserved originals.

### 4.4 Document generation

- Word documents are produced by **executable scripts** (`python-docx` or the Node `docx` library), not hand-assembled.
- Typography precedent from an accepted proposal: Times New Roman body, dark navy H1 headings, charcoal H2 subheadings, black body text.
- Precision layout precedent (3-page proposal): margins 0.65 in top, 0.6 in bottom, 0.7 in sides; body 10.5 pt with line spacing 220 and after-spacing 40; references 9.5 pt with line spacing 220 and after-spacing 20; compressed heading spacing.
- **Page count is a hard requirement, not a target.** Expect to iterate on spacing and font size until it lands exactly.
- Standard proposal section set: Project Overview; Introduction, Problem Statement and Background; Research Questions and Hypothesis; Literature Review; Methodology (with numbered subsections); Clinical Impact and Future Applications; References (5 to 7 real, relevant, verifiable academic citations).
- Citation style for academic work: APA. You have run full APA citation audits and citation gap audits on your own paper.

---

## 5. Established working habits

- You work on multiple parallel projects and expect the assistant to keep them straight by name.
- You iterate. The CV went through repeated passes to a final 3-page Word document. Proposals get restructured after critique.
- You ask for architecture critiques and act on them. Precedent: a six-node fixed-graph GNN design was found to collapse to a linear map, and you took that finding into a restructuring decision rather than defending the original.
- You prepare heavily before meetings, in writing. Meeting-prep documents exist for the Apkarian Lab, CFAAR, and Bhargava Systems.
- You verify claims independently and correct the record when something is filed wrong.
- You separate "what I want built" from "what I want evaluated." The council skill exists specifically for the second.

---

## 6. Active and recent projects

### 6.1 CTPR / Apkarian Lab: drug-receptor database and GNN pipeline
Northwestern, Center for Translational Pain Research, under Dr. A. Vania Apkarian, working alongside Dr. Lejian Huang.

- Building a large-scale drug-receptor binding database plus a GNN-based drug discovery pipeline, aimed at chronic pain drug recommendation.
- Programmatic data pulls from ChEMBL, GtoPdb, PDSP (Ki database), DrugCentral, and DGIdb.
- Primary data gap: functional directionality (agonist versus antagonist annotations).
- Architecture status: the original six-node fixed-graph design collapses to a linear map and needs restructuring. Two options under consideration are a learned-crosstalk GNN on the drug-propagation side, or a bipartite drug-receptor R-GCN with typed edges.
- Runs on Quest (allocation b1090, CPU-only) via VS Code Remote-SSH.
- Public-facing name for this work: "Drug-to-Receptor Space Mapping Tool for Pharmacological Research" (no fMRI mention in that framing).

### 6.2 Chronic back pain opioid R analysis
- R pipeline over a chronic back pain dataset with three groups: healthy, CBP-O (no opioids), CBP+O (long-term opioids).
- Five source spreadsheets: behavior, receptor AI values, PCA results, medication list, MQS. All merged into one master dataset.
- Graph 1 (formerly graph 3): disease axis versus clinical outcome, faceted by opioid subgroup, produced with NRS and with PC1, PC2, PC3 on the Y axis. Former graphs 1 and 2 were dropped.
- Characteristics table modeled on the source paper's Supplementary Table 1, extended with MME, ROE, DOU, PC1, PC2, PC3, MQS, and the count of patients on antidepressants.

### 6.3 B-ALL PAX5/NOTCH1 biomarker meta-analysis
Your foundational independent research project, mentored by Mr. Schramm and Mr. Donahue.

- GEO datasets GSE11877, GSE19475, GSE11504 on platform GPL570, analyzed in JMP Student Edition.
- ANOVA results: PAX5 F(2,302) = 1190.01, eta-squared = .887; NOTCH1 F(2,302) = 2052.35, eta-squared = .931.
- SVM classifier: approximately 97.1% validation accuracy, AUC = 1.000.
- Submitted to IJAS. Gold Award at regionals, advanced to state.
- Full APA citation audit and citation gap audit completed. Publication-quality figures in NIH/journal style.

### 6.4 CFAAR pediatric allergy risk predictor
Northwestern, under Dr. Ruchi Gupta and Dr. Christopher Warren.

- Multimodal pediatric allergy risk predictor.
- Portfolio of computational biology proposals developed for the center: conditional diffusion models, GNNs for allergic march prediction, Bayesian neural networks for personalized oral immunotherapy dose curves, RNA velocity for tolerance prediction, spatial transcriptomics.
- Follow-up sent proposing collaboration with the biostatistics team and the smart inhaler project lead.

### 6.5 RenalSense (Conrad Challenge, Health Sciences)
Teammates: Arnav Nanda, Ishan Suresh Kumar.

- Low-cost flexible chest patch combining single-lead ECG and bioimpedance (MAX30001-class chip) with ML to infer potassium trend and fluid status in dialysis patients. Targets hyperkalemia and fluid overload.
- You led planning through hardware and software architecture, dataset sourcing (PTB-XL, MIMIC-IV-ECG), a GitHub scaffold, and a bill of materials.
- Hardware: MAX30001 breakout, ESP32, MPU-6050, hydrogel electrodes.
- Full 6-month build plan exists.
- Outreach to Professor John A. Rogers at Northwestern, arranged after nephrologist Dr. Lorenzo Gallon recommended the connection.

### 6.6 Medical device concept portfolio
Eight concepts being standardized and vetted for novelty and feasibility, with the constraint that a high school team must be able to execute them.

- **NeuroEcho** (actively pursued as "Team NeuroEcho" with Arnav Nanda): closed-loop sEMG plus edge-AI plus TENS wearable for phantom limb pain in amputees. Concept stage, nothing built. Being pitched to scientists and researchers for feasibility feedback, advisors, IRB guidance, and access to subjects, datasets, and lab hardware.
- **NeuroFlex** (actively pursued as "Team NeuroFlex" with Arnav Nanda): wearable EMG-triggered TENS/NMES band for upper-limb spasticity in stroke survivors. Slide deck presented to academic researchers, clinical neurologists, and PIs for critique. Concept stage.
- **SpastiCast:** iOS/LiDAR scanning to auto-generate 3D-printable custom WHFO splints for stroke and CP spasticity, pitched B2B to OT clinics. (Note: dropped from LinkedIn as of Aug 2026.)
- **AuriCalm:** closed-loop taVNS earbud.
- Levodopa sweat-sensing dose-window patch.
- Electrical impedance tomography concept.
- Muscle NIRS closed-loop concept.
- Forearm A-mode ultrasound concept.
- ALS muscle-tracking concept.

Briefs include technical architecture, reimbursement and business strategy, and pitch demo plans.

### 6.7 NOR behavioral analysis deck
- Novel Object Recognition results deck co-authored with Krish, prepared for Marivi. Compound PP01 versus Vehicle.
- Compares 1-week and 1-month post-treatment, 300 s and 600 s scoring windows. Measures: discrimination index (d2), recognition index (%), total exploration time.
- Sample sizes: Vehicle n = 10 at both timepoints; PP01 n = 10 at 1 week, n = 9 at 1 month (carried from the companion social recognition dataset, pending verification).
- Structural decision made: all Recognition Index slides moved to an appendix because RI is mathematically identical to d2 (RI = 50 x (d2 + 1)), confirmed numerically against plotted values.
- Open items flagged for confirmation before circulation: a "(10m)" duration mislabel on 300-second panels, undefined error bar type in captions, final sample sizes against the NOR animal list, and object identity during habituation.
- Inferential analyses still needed: two-way ANOVA with within-animal structure, one-sample tests against zero, exploration time as covariate, sex as a factor, multiple comparison correction.

### 6.8 Hydrogel wound dressing proposal
- 16-formulation chitosan-alginate hydrogel wound dressing proposal, including BSL-2 safety tables and the full formulation matrix.
- Targeting lab access at Lewis University via the Keleher Research Group (Dr. Jason Rago, graduate researcher Katey Sheets).
- Part of a broader ISEF-level idea exploration covering biomaterials, synovial fluid rheology, and cardiac FEA.

### 6.9 Research outreach system
- Bulk cold outreach to 20+ professors across Northwestern, UChicago, UIC, Loyola, Rush, and DePaul. Earlier rounds also covered Lewis University, North Central College, Benedictine University, and NIU.
- Roughly 15 cold outreach emails drafted and filed as Gmail drafts, with the IJAS meta-analysis paper and abstract attached.
- Labs contacted include: Kiskinis Lab (ALS/iPSC), Cheng Lab (glioblastoma), GMCF at Rush (Dr. Stefan Green), Loyola Stritch faculty, Wong at Northwestern, Dr. Kavouras at Lewis, Dr. González Aparicio and Dr. Visick at North Central.

**Cold email generation rules (standing):**
1. Combine multiple selected project ideas into one pitch.
2. Wet-lab projects take priority in framing over computational ones unless specified otherwise.
3. Search for and include the professor's email address in the final draft.
4. Indirectly hint that the internship is for learning and experience rather than payment, without stating it explicitly.
5. Structure as lab help first, then personal project ideas, with equal emphasis unless specified.
6. All project ideas must tie back to medicine or human health.
7. Computational ideas must be framed as high-impact and transformative, never incremental.
8. **Never propose biomarker discovery as a computational idea.** Aim instead for treatment prediction engines, drug resistance modeling, metabolic pathway simulation, or clinical decision systems.

### 6.10 Bhargava Systems Research proposals
Met with Yash Bhargava and developed three proposals:
- Hybrid GNN-Transformer for polymer-peptide biocompatibility prediction.
- Pediatric off-target toxicity filter for high-throughput virtual screening.
- MeshGraphNet cardiac FEA surrogate.

### 6.11 Pediatric off-target toxicity GNN (proposal detail)
- Secondary screening ML pipeline filtering HTVS hits for pediatric-specific off-target toxicity. Targets NOTCH1, PAX5, CDK4/6, BCL-2 family, leukemia-associated receptors.
- HTVS input scale 10^4 to 10^5 hits (AutoDock Vina or Glide).
- Five multi-task toxicity endpoints: cardiotoxicity, hepatotoxicity, nephrotoxicity, developmental/growth toxicity, neurotoxicity.
- Datasets: Tox21, ToxCast, ChEMBL, DILIrank, hERG-Central, FAERS pediatric subset, with a placeholder for internal lab data.
- Architecture: D-MPNN with ontogeny-conditioned readout as primary; R-GCN heterogeneous knowledge graph as backup/ablation with documented tradeoffs.
- Validation: retrospective against doxorubicin, vincristine, asparaginase, methotrexate; prospective zebrafish embryo wet-lab.

### 6.12 Cardiac FEA surrogate (proposal detail)
- MeshGraphNet trained to predict Cauchy stress tensor redistribution in patient-specific ventricular meshes after simulated surgical incisions (ventriculotomy, septal myectomy, ablation lesion sets).
- Holzapfel-Ogden constitutive model. FEBio as primary solver, SimVascular conditional on licensing.
- Data: public CMR datasets (ACDC, UK Biobank) plus potential institutional data.
- One-year baseline scope structured for multi-year extension, MeshGraphNets in Year 1 with equivariant extensions deferred to Year 2.
- "Topology token" conditioning was flagged as a proposed contribution rather than an established technique.

### 6.13 Earlier computational work
- Cardiovascular biomarker ML pipeline on the Cleveland Heart Disease Dataset (Python/Jupyter).
- SVM classifier built in JMP.

---

## 7. People directory

| Person | Relationship |
|---|---|
| Dr. A. Vania Apkarian | PI, Apkarian Lab / CTPR, Northwestern. Proposals submitted: snRNA-seq OIH stratification, simultaneous fMRI + fiber photometry, chemogenetic MOR+ NAc silencing, pre-op brain scan opioid prediction model. |
| Dr. Lejian Huang | Apkarian Lab / CTPR. Involved in meeting prep and follow-ups. |
| Dr. Ruchi Gupta | Director, CFAAR, Northwestern. |
| Dr. Christopher Warren | CFAAR researcher. You work under him on the allergy predictor. |
| Mr. Schramm | Research mentor, Waubonsie Valley. Co-mentored the B-ALL meta-analysis. |
| Mr. Donahue | Research mentor, Waubonsie Valley. Co-mentored the B-ALL meta-analysis. Reference for the YAC Science Olympiad application. |
| Arnav Nanda | Teammate on RenalSense, Team NeuroEcho, Team NeuroFlex. |
| Ishan Suresh Kumar | Teammate on RenalSense. |
| Prof. John A. Rogers | Northwestern. Contacted for RenalSense. |
| Dr. Lorenzo Gallon | Nephrologist. Recommended the Rogers connection. |
| Dr. Jason Rago | Lewis University, Keleher Research Group. Hydrogel project contact. |
| Katey Sheets | Graduate researcher, Lewis University, KRG. Hydrogel project contact. |
| Dr. Stefan Green | GMCF (Genomics and Microbiome Core Facility), Rush. |
| Yash Bhargava | Bhargava Systems Research Inc. Three proposals developed. |
| Krish | Co-author on the NOR deck. |
| Marivi | Recipient of the NOR deck. |
| Hasan Al-Khalidi | Prepared the social recognition deck used as your house style reference. |

---

## 8. Domain shorthand glossary

So a new assistant does not have to ask:

- **CTPR** = Center for Translational Pain Research (Apkarian Lab, Northwestern)
- **CFAAR** = Center for Food Allergy & Asthma Research (Northwestern)
- **CBP-O / CBP+O** = chronic back pain without / with long-term opioids
- **MQS** = Medication Quantification Scale; **MME** = morphine milligram equivalents; **ROE** = rate of opioid escalation (verify expansion); **DOU** = duration of opioid use (verify expansion)
- **NRS** = numeric rating scale (pain)
- **NOR** = Novel Object Recognition; **d2** = discrimination index; **RI** = recognition index
- **B-ALL** = B-cell acute lymphoblastic leukemia
- **HTVS** = high-throughput virtual screening
- **OIT** = oral immunotherapy; **OIH** = opioid-induced hyperalgesia
- **KRG** = Keleher Research Group (Lewis University)
- **IJAS** = Illinois Junior Academy of Science; **BPA** = Business Professionals of America; **SNHS** = Science National Honor Society
- **WHFO** = wrist-hand-finger orthosis
- **taVNS** = transcutaneous auricular vagus nerve stimulation
- **Quest** = Northwestern's HPC cluster

---

## 9. Academic and extracurricular record

**Coursework (recorded):** AP Biology, AP Language and Composition, honors precalculus/calculus, Anatomy & Physiology, Spanish 4. (verify: this list is from an earlier term and is likely out of date given the College of DuPage dual enrollment.)

**Activities and honors:**
- Varsity tennis.
- Science Olympiad: state medalist, Paula Mueller Leadership Award. Applied to the Waubonsie exec board and the Illinois Science Olympiad Youth Advisory Committee state board.
- BPA: national finalist in Health Research Presentation and Computer Security. Chapter board application built around a platform called "The Blueprint" (quarterly event series, mentorship pairing, career panels).
- SNHS: board member since August 2026. Service Officer application platform was "Lab Explorers" (youth outreach, career-forward programming, Fermilab/Argonne/hospital partnerships).
- IJAS: State Competition Gold Award. Also applied for the IJAS research mentor role.
- Key Club, AMP Club, Speech Team.
- Lead volunteer at Arbor Terrace Senior Living. Nursing home volunteer leadership is a recurring evidence point in your applications, alongside the B-ALL research.
- Ran a high school Science Club, planning experiments and design-build competitions.
- Science Olympiad Protein Modeling: protein 2LV8 (PDB: OR16), Jmol command sets, physical Mini-Toober build, Model Description Legend, NAD+ binding pocket as the creative addition. Competed at state.
- **Not** involved with the District Student Academy Board (this was recorded incorrectly once and corrected).

**Applications submitted:** UIUC SpHERES Young Scholars, NURPH (Northwestern materials science summer program), YAC Science Olympiad executive board.

**Academic projects:** AP Biology transpiration lab (CER with class-wide data), trophic cascades worksheet (HHMI Click & Learn), AP Bio poster on exa-cel/Casgevy (sickle cell CRISPR trial, Frangoul et al. 2024 NEJM), Anatomy & Physiology semester review, AP Lang synthesis essay on animal research ethics, AP Bio population ecology FRQ tool, thermos insulation design challenge, Spanish Socratic seminar prep and scripts.

**Professional materials:** CV built iteratively to a final 3-page Word document covering three research placements (Apkarian Lab CTPR, CFAAR, RenalSense/Conrad). Teacher bragsheet built for recommendation letters. LinkedIn headline: "Research Intern @ NU FSM | Biomedical/Clinical Research & Device Design | Student @ Waubonsie Valley | Aspiring Physician-Scientist."

**LinkedIn-specific rules (in case you automate it again):** omit dual enrollment references entirely; list broad skills (Data Science, Neural Networks, Deep Learning) rather than specific tools (ggplot2, lavaan, XGBoost, Git, GitHub, VS Code, Miniforge); IJAS honor reads "Gold Award," not "Gold Medalist"; SNHS goes under Organizations/Leadership, not Volunteering. Known platform quirk: LinkedIn auto-rewrites a custom headline when a new current-position experience entry is saved.

---

## 10. Interests outside work

Sneakers, fragrances, PC gaming (Legion 5 Pro), fantasy football, photography (iPhone ProRAW to Lightroom workflow), graphic design.

---

## 11. Trajectory: past, present, future

**Past:** Started with independent high school research, the B-ALL PAX5/NOTCH1 meta-analysis, mentored in-house at Waubonsie Valley. Built early computational chops on the Cleveland Heart Disease dataset and JMP-based classification. Ran a systematic cold outreach campaign across Chicago-area universities to convert that work into lab access.

**Present (as of Aug 2026):** Three concurrent research affiliations (CTPR/Apkarian, CFAAR, plus the Mind4Youth editorial role). Two active device teams with Arnav Nanda. Real HPC access. A GNN pipeline mid-restructure. Eight device concepts under novelty and feasibility review.

**Future:** Academic cardiothoracic surgeon combining clinical care, research, and teaching. The near-term pipeline is publication, competition placement (IJAS, Conrad, BPA, Science Olympiad), and converting device concepts from paper into built prototypes with advisor and IRB support.

**One side thread:** you explored entry-level IT and data center careers (Chicago corridor, Elk Grove Village and Aurora) and air traffic control, focused on paths enterable with under a year of training, with willingness to pay for that training. **Confirmed: this research was for a family member, not a change of direction for you.** Do not treat it as a signal about your own goals.

---

## 12. Quick behavioral checklist for the new assistant

Before sending any response:
- [ ] No em dashes anywhere.
- [ ] No "furthermore," "moreover," "in conclusion," "it is important to note."
- [ ] Every factual claim is verified or explicitly marked uncertain.
- [ ] No invented citations, statistics, or API details.
- [ ] Sources cited inline or as a short end note.
- [ ] Nothing existing was deleted or overwritten without asking first.
- [ ] For large deliverables: clarifying questions asked and answered before generating.
- [ ] For generated files: validated programmatically (page count, prohibited phrases, render check) before delivery.
- [ ] Slides delivered as editable .pptx, not PDF.
- [ ] Full name written as "Kushagra Verma" in formal output.

---

## 13. Known gaps (fill these in yourself)

I do not have recorded answers for these, and I am not going to guess:

1. Line length, indentation, and formatter preferences (Black? Ruff? Prettier? tidyverse style?).
2. Type-hint policy in Python, and whether you want docstrings and in what style.
3. Comment density preference.
4. Testing habits (pytest, testthat, or none).
5. Preferred project directory structure and file naming conventions.
6. Whether you prefer notebooks or scripts for exploratory work.
7. Branch naming and commit message conventions.
8. Preferred error handling style.
9. Whether R work uses tidyverse or base R idioms by default.
10. Environment management preference (conda/mamba vs venv vs uv).

If you answer those ten once, this document becomes considerably more useful to any assistant reading it.

---

## Appendix A: Source note

Everything above comes from stored memory files and past conversation records in Claude, current as of August 28, 2026. No external sources were used. No details were inferred beyond the lines explicitly marked `(inferred)`.

---

## Appendix B: Condensed version for a system prompt

> I'm Kushagra Verma (Kush), a high school student in Naperville, IL, dual-enrolled at College of DuPage, working as a research intern at Northwestern Feinberg (Center for Translational Pain Research / Apkarian Lab) and with CFAAR. Long-term goal: academic cardiothoracic surgeon. I work in Python and R on computational biology and ML pipelines (GNNs, drug-receptor modeling, clinical data analysis), plus biomedical device concepts.
>
> Rules for working with me:
> 1. Never use em dashes. Ever. Use parentheses or separate sentences.
> 2. Never use "furthermore," "moreover," "in conclusion," "it is important to note."
> 3. Never fabricate data, statistics, citations, or API details. If uncertain, say so explicitly. Verify library details against current docs rather than memory.
> 4. Cite sources non-disruptively (inline links or a short end note).
> 5. Be candid, not optimistic. Tell me where things break down.
> 6. Only add, never delete or overwrite my existing content without asking first. Bring proposed changes to me before executing.
> 7. For anything large: research first, ask me clarifying questions, then stop and wait. Do not generate a full draft on assumptions.
> 8. Validate generated files before delivering (page count, prohibited phrases, render check). Slides as editable .pptx, never PDF.
> 9. Use "Kushagra Verma" as my full name in formal documents.
> 10. Figures: NIH/journal style, white background, Arial, ColorBrewer, y-axis starts at zero.
> 11. For data analysis: merge sources into one master dataset that all figures read from, show how every plotted value is calculated, and display group means on every graph.

---
---

# ADDENDUM: Second-Pass Findings
**Added:** August 28, 2026

The first version of this document was built from stored memory files. This addendum comes from a second sweep across the full conversation history, run as five parallel passes: coding and tooling, recent activity, communication and writing, research methodology, and personal or biographical detail. Everything here is new material the first pass missed, plus one correction.

---

## A. Correction to the main document

**Section 11, side thread:** confirmed by you. The entry-level IT, data center, and air traffic control research was for a **family member**, not for you. The main text has been updated.

---

## B. How you actually work with code (important, and missing from v1)

This changes how an assistant should treat you, so put it near the top of anything you paste into a new tool.

- **You self-identify as a "vibe coder."** Your stated workflow: describe detailed architecture to the model and have it generate the code, rather than writing the code yourself. You are explicit that this is the preferred mode, not a limitation you are apologizing for.
- **Self-assessed coding level: beginner.** Stated directly when scoping the RenalSense team (one intermediate coder, you as a beginner, one member with no coding experience, two members with ML experience). The team was framed as "AI-Assisted Developers" relying on AI tools to write, debug, and compile.
- **What this means for an assistant:** do not hand back a fragment and assume you will wire it in. Give complete, runnable files. Explain what a change does and why, not just the diff. Do not assume familiarity with build tooling, environment quirks, or language idioms. At the same time, **do not talk down about the domain science**, where your knowledge is well above the coding level. The gap between those two is the thing to calibrate to.
- **Your architectural judgment is strong and should be engaged seriously.** You caught that NeuroEcho and NeuroFlex shared nearly identical hardware stacks, electrode approaches, artifact problems, IRB gates, and surrogate-validity limits, and that a feasibility table distinguishing them as "Medium" versus "Medium-high" was therefore unjustified. That critique was correct and the table was fixed. Expect this and welcome it.
- **You ask to be criticized, explicitly.** Recorded instruction on the GNN work: do not simply agree, identify gaps, be critical. Recorded instruction on device feasibility: address both positives and negatives, not just strengths.

---

## C. Your full AI tooling stack (as of August 28, 2026)

This was completely absent from v1 and is the most directly useful section for the migration.

**Current daily setup:**

| Surface | Tool | Backend |
|---|---|---|
| Terminal, primary | `claude` (Claude Code CLI) | Anthropic |
| Terminal, Google | `agy` (Antigravity CLI) | Your existing Gemini subscription |
| Terminal, free/experimental | `opencode` | NVIDIA NIM (build.nvidia.com) |
| VS Code sidebar | Claude extension + Antigravity extension | both |
| Positron sidebar | Claude extension + Kilo Code | Claude + NVIDIA |
| Positron terminal | `agy` | Antigravity |

**Editors:** VS Code and **Positron** (the Posit R/Quarto IDE). Positron was missing from v1 entirely. Its value to you is the Console, Variables, and Data Explorer panes, which is why you generally do not want an agent panel eating horizontal space there.

**You work in Quarto.** Also missing from v1. Relevant because open models hosted on NVIDIA are noticeably weaker at R and Quarto than Claude is, since R is underrepresented in most training mixes. The working split settled on: NVIDIA for Python, bulk refactors, and throwaway experiments; Claude for the R analysis work.

**Accounts and keys in play:** Anthropic (Claude Pro), a Google/Gemini subscription, an NVIDIA Developer Program account for build.nvidia.com (free tier, key prefix `nvapi-`, bounded credits and a rate limit, agentic loops burn 30 to 80 calls per task), and Codex was used during RenalSense.

**Ecosystem facts you already know, so an assistant should not re-explain them:** Gemini CLI was deprecated in June 2026 and replaced by Antigravity. Roo Code shut down May 15, 2026 and its repo is archived. Kilo Code rebuilt its CLI on the OpenCode core. Microsoft Marketplace terms restrict `.vsix` use in third-party editors, which is why Cursor, Windsurf, and Antigravity default to Open VSX, and why Positron extension availability is an ongoing friction point.

---

## D. Antigravity-specific setup notes (carry these over)

Direct analogues to your existing Claude Code configuration:

- **`AGENTS.md` at project root.** This is the direct equivalent of your `CLAUDE.md`. Anything in it is prepended to every prompt processed inside that directory. This document is what should go there.
- **`.agents/skills/*.md`** for reusable slash commands. A file at `.agents/skills/lint.md` becomes `/lint` in the TUI. Skills can carry their own prompts and allowed-tools lists.
- **Global skills** live at `~/.gemini/antigravity-cli/skills/`. This is where a rebuilt `council` skill would go.
- **`mcp_config.json`** carries over MCP servers from your Claude Code setup.
- **`agy plugin import gemini`** migrates old Gemini CLI extensions.
- **Two commands with no Claude Code equivalent:** `/goal` keeps a task running to completion without stopping for input, and `/schedule` runs an instruction later, one-time or recurring. Browser use is opt-in behind `/browser`.
- **Watch quota when using subagents.** Directly relevant given the council pattern spawns four in parallel.

**Rebuilding the council in Antigravity:** the four subagent definitions and the orchestrator SKILL.md are reproduced in structure in Section 3.4 of the main document. The read-only tool restriction (no Write, Edit, or Bash) is the part worth preserving deliberately, since it is what makes the council safe to run on anything.

**One caveat you were given and should keep in mind:** this ecosystem churned hard through 2026. The durable filter is not GitHub stars, it is license, last commit date, and how the project makes money. Because you wire these to your own keys rather than a bundled subscription, switching harnesses costs a config file, not a workflow.

---

## E. Communication rules that v1 missed

These are as binding as the ones in Section 2.

### E.1 Punctuation in your own voice

When editing or proofreading **your** writing, avoid **colons, em dashes, and semicolons**. You stated plainly that those do not reflect how you write. This is broader than the global em dash ban and applies to anything written in your voice: application essays, About sections, personal statements, emails.

### E.2 Proofreading protocol

- **Always return a numbered change log.** Every correction listed with the original text and the revised text. No silent edits, ever.
- **You draft, the assistant proofreads.** Your strong preference is to write your own version and have it copyedited, rather than receiving a generated draft. When you do want a draft, you provide raw thoughts first and have those shaped, rather than receiving something unsolicited.
- **Light copyedits only by default:** spelling, grammar, capitalization, minor phrasing. Not rewrites.
- **Preserve your voice.** You have rejected generic phrasing repeatedly (for example, an opening about being "most confident" in a quality was rejected as unnatural).

### E.3 Inline over files

**Deliver responses inline in chat rather than as separate files unless a file is actually necessary.** This was stated explicitly during application work. Files are for real artifacts (decks, documents, scripts), not for text you are going to read once and paste somewhere.

### E.4 Register switching

For non-technical or casual questions, you want **short and plain**. Recorded feedback on an over-long, over-formal answer: keep it simple, casual, no big or weird words, one to two sentences per point rather than paragraphs. Clear structure, no academic tone. The technical rigor in Section 2 applies to research and code work, not to everything.

### E.5 Cross-references over repetition

When content overlaps across sections of a multi-part application or document, flag the overlap and handle it with a cross-reference rather than repeating yourself. You made structural calls like keeping a Vice President statement short and pointing reviewers back to the President paragraph for shared qualities.

### E.6 Audience realism

You have good instincts about what is common in a given applicant pool and adjust emphasis accordingly. Example: you noted that research projects are common among IJAS applicants and should therefore not be overemphasized in an IJAS application. An assistant should not assume your strongest credential is always the one to lead with.

### E.7 Explicit gates

You set your own go/no-go signals and expect them to be honored. During the LinkedIn session you told the assistant to wait for the phrase "I am in" before touching the live profile, and the assistant held. Watch for gates like this and do not proceed past them.

---

## F. Research methodology preferences (new)

### F.1 Adversarial literature search

This is your standout methodological rule and it was missing entirely from v1. When evaluating a concept, do not run confirming searches. Run **second-order searches that actively seek disconfirming evidence**: null results, non-significant trials, author overlap between supposedly independent papers, and commercial incumbents already doing the thing.

Findings this protocol has actually produced for you:
- The Tai 2019 and Guo 2026 levodopa papers share a senior author.
- The PNAS 2026 Saha paper closed most remaining gaps in levodopa sweat sensing.
- The "world's first closed-loop" taVNS claim for AuriCalm was contradicted by at least four prior systems.
- Atanackov 2025 RCT found no significant RMSSD change at the cymba conchae.
- Lendaro 2025 multicenter RCT found phantom motor execution non-superior to imagery.
- ActivArmor already commercializes the SpastiCast pipeline, including the identical reimbursement argument.
- Monte-Silva 2019 gives SMD 0.63 at impairment level in chronic stroke but no activity-level effect; a 2022 meta-analysis found SMD 0.14, not significant versus non-sEMG interventions.
- The commonly cited 80% phantom limb pain prevalence figure does not hold. The pooled meta-analytic estimate is 64% (Limakatso et al. 2020), with a lifetime range of 76 to 87% and one-year incidence of 82% (Stankevicius et al. 2021).

You explicitly said the value in these responses came from the second-order searches, not the initial confirming ones.

### F.2 Claim discipline

- **Separate engineering claims from clinical claims.** When the literature does not support a therapeutic mechanism claim, reposition the contribution as an engineering claim (closed-loop at wearable scale and cost) and put the contested evidence on its own dedicated slide rather than burying it.
- **State what you lack.** Your NeuroEcho deck explicitly listed what the team does not have: IRB standing, clinical access, credentialed oversight. Keep that habit.
- **Caveat borrowed numbers.** Example flagged: the 1 to 2 percent accuracy loss figure for post-training quantization is not sEMG-specific and must be stated with that caveat.
- **Replace figures that will not survive scrutiny**, even when they are widely cited, and label the replacement with its source.

### F.3 Standardized device brief format (eleven sections)

Your house format for concept briefs, used across all eight device concepts:

1. Concept statement
2. Clinical problem
3. Mechanism
4. Proposed solution
5. Technical architecture and BOM
6. Prior art
7. Market structure
8. Regulatory classification
9. Demonstration requirements
10. Risks
11. Evidence assessment

Rules for these briefs, which you corrected the assistant into: **no opinions, no verdict sections, no second-person address, no references to you or your team, no cross-concept links.** Purely objective and informative. Each brief must be self-contained, because they are fed to the council for independent review. Inline author-year citations with per-concept reference lists. A front matter section explains document scope and method.

**Council scaling note:** for multi-concept review, one agent per concept works better than one agent processing all concepts at once, paired with an adversarial search protocol.

### F.4 Regulatory and domain precision

Corrections you have absorbed and now expect to be handled correctly: the distinction between FDA-regulated stimulation devices and student-built circuits, correct use of the term "SaMD," CPT/HCPCS billing prerequisites, and realistic reimbursement figures.

---

## G. Projects and artifacts missing from v1

- **The master analysis script:** a 4,394-line R script in RStudio using ggplot2 and lavaan, for statistical modeling and figure generation across clinical datasets. This is the concrete scale of the CBP opioid pipeline described in Section 6.2. Any assistant working on it should know it is a large single script, not a small one.
- **NeuroFlex deck build (Aug 2026):** 17 slides via python-pptx, with speaker notes and citations on every slide, an evidence slide separating supported from contested claims, and discussion questions for PI reviewers. Built with a render-and-inspect QA loop (build, validate, convert to PDF via LibreOffice headless, rasterize with pdftoppm at 96 to 100 dpi, view individual slide JPEGs).
- **NeuroEcho deck:** 15 slides, house style extracted from three of your existing lab meeting decks. Palette navy `0B2545` and teal `0EA5A5`, Cambria headings, Calibri body, numbered eyebrow labels, takeaway bands. Risk slide ordered by project-killing severity.
- **NeuroFlex 3D concept render:** a self-contained Three.js HTML file showing an assembled and exploded five-layer model (forearm reference, electrodes and stimulation pads, conductive liner, elastic strap, control module housing). Anatomy drove layout: sEMG bipolar pairs volar over the flexor compartment, NMES pads dorsal over the extensors in two axial rows for field steering, housing radial and clear of pads. Geometry was verified numerically in Node rather than by screenshot.
- **NeuroFlex ML rationale (settled):** LDA and linear SVM as baselines because LDA is the field-standard baseline in myoelectric pattern recognition with a closed-form solution and minimal compute. The held-out-session gate is the load-bearing constraint, since within-session accuracy is misleading. TCN as the escalation target for deterministic latency via dilated causal convolutions. CMSIS-NN kernels are the mechanism behind speedups beyond the 4x from int8 quantization. Two open concerns: the literature is dominated by prosthesis control in able-bodied and amputee populations rather than post-stroke forearm EMG, and per-donning calibration is a genuine usability problem for a stroke survivor working one-handed on the affected arm. Version one should trigger and then pause classification during the stimulation burst rather than recording concurrently. **PhysioMio** was identified as a strong direct dataset fit.
- **MacBook Air vertical dock (personal 3D printing project):** a five-part parametric assembly designed in **build123d** and exported to STEP and STL, intended for further editing in **Fusion 360**. Auto-mates with a UGREEN Revodok Pro 106 USB-C hub. Design philosophy: the laptop bottoms out on a printed slot floor rather than loading the USB-C port, with a floating carrier giving 2.5 mm vertical travel and 1.5 mm lateral float per side. Four measurements are still placeholders pending your calipers: `PORT_OFFSET` (60 mm placeholder), `BOOT_W`/`BOOT_D`, `SHELL_PROT`, and `INSERT_D`. Tightest constraint is a 0.5 mm seating margin.
- **Additional research proposals not listed in v1:** a conditional diffusion model for pediatric airway remodeling, and a spatial transcriptomics CNN related to virtual staining.
- **A medical prefix and suffix reference document** built in Word.
- **IJAS 100th anniversary theme work:** prior year's theme was "Driven by Discovery: Research is Key," which sets the format you want (concise, dual-clause, inspirational).
- **Contact lookup precedent:** Prasad Shirvalkar at UCSF, `Prasad.Shirvalkar@ucsf.edu`, sourced from a ClinicalTrials.gov protocol document. You ask for sourced results on contact lookups, not bare answers.

---

## H. Background detail missing from v1

**Certifications:** CPR. Teen Mental Health First Aid (tMHFA). Certiport Cybersecurity Verified.

**Additional activities:** STEMSTERS (out-of-school STEM outreach where you teach younger students). American Red Cross. Shadowing at Athletico Physical Therapy.

**Additional honors on record:** Science Olympiad State, Gene Regulation, 2nd place. Science Olympiad State awards in Health Science and Experimental Design. Scullen Star Leadership Award. Top Science Student Award. Student of the Semester. Paula Mueller Leadership Award (already in v1).

**One discrepancy to resolve yourself:** your LinkedIn lists BPA Nationals **Finalist** for Computer Security but BPA Nationals **Qualifier** for Database Applications, while your written brief said Finalist for both. Confirm which is correct before it goes on an application.

**Education history:** Waubonsie Valley High School, 2024 to 2028. Scullen Middle School before that. Dual enrollment at College of DuPage (omitted from LinkedIn by your own instruction).

**Currently in calculus.** Recent work covers limits, one-sided limits with absolute values, limit laws, difference quotients, infinite limits and vertical asymptotes, limits at infinity, derivatives, and tangent lines through external points. Your pattern is to attempt problems independently first and then ask for verification, and to ask targeted follow-ups on single intermediate steps rather than whole problems.

**The nursing home story** is your recurring leadership anecdote and it is worth keeping straight: at Arbor Terrace you noticed residents were being sidelined during art activities, advocated to a supervisor for a new approach (hand-over-hand guidance, resident-directed conversation, adapted materials), and the change increased resident engagement and pride. It has been used across the IJAS board, Science Olympiad board, and SNHS applications.

**Personal:** you celebrate Rakshabandhan and have cousins you write to for it.

**Known critiques you have of IJAS State Exposition** (used in board applications): judging conflicts of interest, unpredictable board-standing times, limited food infrastructure, car dependency, and weak connections to the broader STEM ecosystem. Proposed fixes included staggered judging, conflict-of-interest reform, university and vendor partnerships, and partnerships with BioGenius, the Chicago Youth Innovation Conference, and MATTER.

---

## I. Accounts and identifiers

- Northwestern email: `kushagra.verma@northwestern.edu`
- LinkedIn: `linkedin.com/in/kushagra-verma01`
- Notion workspace ID: `d6867e37-eeae-42cd-87fa-265f511a1575`. Connection can be verified with a self-fetch (`id: self`).
- Quest NetID: recorded as `okt1985` (still flagged, verify before use)
- Quest allocation: `b1090`, node `qnode3027`

---

## J. People to add to the directory

| Person | Relationship |
|---|---|
| Dr. Rubens | Mini Medical School. You applied to be a Teaching Assistant, submitting a curriculum and station-improvement proposal. |
| Prasad Shirvalkar | UCSF. Contact looked up for outreach, pain and neuromodulation research. |
| Your Science Olympiad coach | Also the YAC sponsor and application reviewer. Worth remembering that these are the same person when writing anything for either. |

---

## K. LinkedIn incident worth remembering

During the Playwright session, saving a new current-position experience entry caused LinkedIn to auto-overwrite your custom headline, destroying the original. The assistant had recorded the old text earlier in the conversation so it was recoverable, and it told you immediately rather than letting you discover it later.

Two lessons that generalize to any automation an assistant runs on your accounts:
1. **Capture the current state before any write, not after.**
2. **Report self-caused damage immediately and plainly**, including the fact that it violated your standing rule, rather than quietly fixing it.

---

## L. Revised gap list

Resolved by this pass: your editors (VS Code and Positron), your language mix including Quarto, your self-assessed coding level, your AI tooling stack, your proofreading protocol, and your inline-versus-file preference.

Still genuinely unknown, and worth answering once:

1. Formatter and linter preferences (Black, Ruff, Prettier, styler, air).
2. Whether R work uses tidyverse or base R idioms, and whether you prefer the native pipe `|>` or `%>%`. (An `AGENTS.md` template drafted for you contained "Prefer tidyverse idioms, use `|>` not `%>%`" but that was a suggested example, not your stated rule. Do not let a future assistant treat it as yours.)
3. Whether you use `renv` for R project environments.
4. Type hints and docstrings in Python.
5. Testing habits (`pytest`, `testthat`, or none).
6. Project directory structure and file naming conventions.
7. Branch naming and commit message conventions.
8. Environment management (conda, mamba, venv, uv).
9. Whether exploratory work goes in notebooks, scripts, or Quarto documents.
10. Preferred error handling style.
