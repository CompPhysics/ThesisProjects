%% University of Oslo thesis document class

%% Identification:

\NeedsTeXFormat{LaTeX2e}[1995/12/01]  % Any version, really!
\ProvidesClass{uiomasterthesis}
  [2023/09/15 v1.09 UiO document class for a master's thesis]

%% Required packages:
\RequirePackage[T1]{fontenc}
\RequirePackage{xifthen}

%% Option declarations:
\DeclareOption{binding}{\def \uiom@bind {24pt}}
\DeclareOption{font=cmr}{\def \uiom@font {\RequirePackage{lmodern}
    \RequirePackage[scaled]{helvet}}}
\DeclareOption{font=garamond}{\def \uiom@font {\RequirePackage[garamond]{mathdesign}
    \RequirePackage[scaled]{helvet}}}
\DeclareOption{font=noto}{\def \uiom@font {\RequirePackage{noto}}}
\DeclareOption{font=times}{\RequirePackage{txfonts}}

% Default values:
\def \uiom@bind {0pt}

% Process options:
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ExecuteOptions{font=cmr}
\ProcessOptions \relax

%% Load base class:
\LoadClass[a4paper,11pt,twoside,openright]{report}

%% Language adaptions:
\newcommand{\papername}{Paper}

%% Margins:
\RequirePackage[vmargin=3cm,bindingoffset=\uiom@bind]{geometry}

%% Text adjustments:
\pretolerance = 3000
\tolerance = 6000   \hbadness = \tolerance
\setlength{\parskip}{0cm plus 1mm}

%% Select a font
\uiom@font

%% Chapter and section headers
%% (mainly to use \sf for these):
\RequirePackage{titlesec}
\titleformat{\part}[display]
  {\centering\sf\bfseries\Huge}{%
    \partname~\thepart}{15pt}{\Huge\sf\bfseries}
\titleformat{\chapter}[display]
  {\sf\huge\bfseries}{\chaptertitlename\ \thechapter}{15pt}{\Huge\sf\bfseries}
\titlespacing*{\chapter} {0pt}{50pt}{40pt}
\titleformat{\section}
  {\sf\Large\bfseries}{\thesection}{1em}{}
\titlespacing*{\section} {0pt}{3.5ex plus 1ex minus .2ex}{2.3ex plus .2ex}
\titleformat{\subsection}
  {\sf\large\bfseries}{\thesubsection}{1em}{}
\titlespacing*{\subsection} {0pt}{3.25ex plus 1ex minus .2ex}{1.5ex plus .2ex}
\titleformat{\subsubsection}
  {\sf\normalsize\bfseries}{\thesubsubsection}{1em}{}
\titlespacing*{\subsubsection}{0pt}{3.25ex plus 1ex minus .2ex}{1.5ex plus .2ex}
\titleformat{\paragraph}[runin]
  {\sf\normalsize\bfseries}{\theparagraph}{1em}{}
\titlespacing*{\paragraph} {0pt}{3.25ex plus 1ex minus .2ex}{1em}
\titleformat{\subparagraph}[runin]
  {\sf\normalsize\bfseries}{\thesubparagraph}{1em}{}
\titlespacing*{\subparagraph} {\parindent}{3.25ex plus 1ex minus .2ex}{1em}

%% Define \uiopaper and modify the numbering
%% in \part-s (only in uiophdthesis).
\newif \ifuiom@in@paper  \uiom@in@paperfalse

%% Insert printed material:
\RequirePackage{pdfpages}
\newcommand{\uioincludepdf}[2][]{{%
  \let \ps@none = \ps@empty
  \def \ps@default {\let\@mkboth\@gobbletwo
    \let\@oddhead\@empty
    \def\@oddfoot{\reset@font\hfil\textsf{\thepage}\hspace*{-2cm}}%
    \let\@evenhead\@empty
    \def\@evenfoot{\reset@font\hspace*{-2cm}\textsf{\thepage}\hfil}}%
  \def \ps@low {\let\@mkboth\@gobbletwo
    \let\@oddhead\@empty
    \def\@oddfoot{\reset@font\hfil\raisebox{-1.25cm}{\textsf{\thepage}}\hspace*{-2cm}}%
    \let\@evenhead\@empty
    \def\@evenfoot{\reset@font\hspace*{-2cm}\raisebox{-1.25cm}{\textsf{\thepage}}\hfil}}%
  \def \uiom@ps {default}
  \setkeys{uiom@keys}{#1}
  \cleardoublepage
  \includepdf[pages={-},pagecommand={\thispagestyle{\uiom@ps}}]{#2}}}
\define@key{uiom@keys}{numbers}{\def \uiom@ps {#1}}

%% Redefine the table of contents and lists of figures and tables:
\RequirePackage{titletoc}
\contentsmargin{0.9cm}
\titlecontents{part}[0.0em]{\addvspace{2.25em plus 1pt}\large\sf\bfseries}{}{}
  {\hfill\contentspage}
\titlecontents{chapter}[2.5em]{\sf}{\contentslabel{2.5em}}{\hspace*{-2.5em}}
  {\titlerule*[1pc]{.}\contentspage}
\titlecontents{section}[5.5em]{\sf}{\contentslabel{3.0em}}{\hspace*{-3.0em}}
  {\titlerule*[1pc]{.}\contentspage}
\titlecontents{subsection}[9.0em]{\sf}{\contentslabel{3.5em}}{\hspace*{-3.5em}}
  {\titlerule*[1pc]{.}\contentspage}
\titlecontents{subsubsection}[13.0em]{\sf}{\contentslabel{4.0em}}{\hspace*{-4.0em}}
  {\titlerule*[1pc]{.}\contentspage}
\titlecontents{paragraph}[17.5em]{\sf}{\contentslabel{4.5em}}{\hspace*{-4.5em}}
  {\titlerule*[1pc]{.}\contentspage}
\titlecontents{subparagraph}[22.5em]{\sf}{\contentslabel{5.0em}}{\hspace*{-5.0em}}
  {\titlerule*[1pc]{.}\contentspage}
\titlecontents{figure}[3.0em]{\sf}{\contentslabel{3.0em}}{\hspace*{-3.0em}}
  {\titlerule*[1pc]{.}\contentspage}
\titlecontents{table}[3.0em]{\sf}{\contentslabel{3.0em}}{\hspace*{-3.0em}}
  {\titlerule*[1pc]{.}\contentspage}

%% Redefine \begin{abstract}...\end{abstract} and
%% define \begin{xabstract}[heading] ... \end{xabstract}
% First, we need a modified quotation environment
% with no indentation of the first line:
\newenvironment{uiom@noindent@quotation}
  {\list{}{%
     \listparindent = 1.5em\relax
     \rightmargin = \leftmargin
     \parsep \z@ \@plus\p@}%
     \item\relax}
  {\endlist}
% Then, we need an environment for formatting the abstract:
\newenvironment{uiom@abstract}[1][\abstractname]
  {\begin{center}\small\textsf{\textbf{#1}}\end{center}\uiom@noindent@quotation}
  {\enduiom@noindent@quotation}             
\renewenvironment{abstract}[1][\abstractname]
  {\cleardoublepage \pagestyle{empty}\vspace*{5\baselineskip}
    \uiom@abstract[#1]}
  {\enduiom@abstract}
% We also define an environment for a second abstract:
\newenvironment{xabstract}[1][\abstractname]
  {\par\vspace{3\baselineskip}\uiom@abstract[#1]}
  {\enduiom@abstract}


%% Define \begin{preface}
\newenvironment{preface}{%
  \chapter*{\prefacename}
  \markboth{\prefacename}{}
  }
  {}

%% Headers and footers:
\RequirePackage{fancyhdr}
\fancyhead{}
\setlength{\headheight}{15pt}
\fancyhead[LE]{\textsf{\small \nouppercase{\leftmark}}}
\fancyhead[RO]{\textsf{\small \nouppercase{\rightmark}}}
\renewcommand{\headrulewidth}{0pt}
\fancyfoot{}
\fancyfoot[LE,RO]{\textsf{\thepage}}
\renewcommand{\footrulewidth}{0pt}
\pagestyle{fancy}
\fancypagestyle{plain}{\fancyhead{}}

%% Figure and table captions
\RequirePackage[labelfont=sf,font=small]{caption}

%% Document parts:
\providecommand{\frontmatter}{\cleardoublepage \pagenumbering{roman}}
\providecommand{\mainmatter}{\cleardoublepage \pagenumbering{arabic}}
\providecommand{\backmatter}{\cleardoublepage}

\endinput
