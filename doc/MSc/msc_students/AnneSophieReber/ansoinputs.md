Hei Anso, nå har jeg gått gjennom det du har skrevet og syns det lover veldig godt.

Du har ei fin framtislling av de relevant3 vitenskapelige spørsmål med
et bra praktisk grunnlag: innhenting av databaseoppføringer,
konstruksjon av deskriptorer, undersøkelse av manglende data,
sammenstilling av litteraturbaserte klasseetiketter og utforskning av
det resulterende kjennetegnsrommet. Utkastet identifiserer også flere
problem, blant annet usikre negative etiketter og forskjeller mellom
emisjonsmekanismer. Dette er gode utgangspunkter for dine studier og
sjølve oppgava.

Manuskriptet viser (så langt jeg klarte å få med meg) foreløpig ikke
en validert metode for å identifisere nye vertsmaterialer for
enkeltfotonemittere. De viktigste hindringene er etikettenes betydning
og uavhengighet, validering på materialer som faktisk er ukjente for
modellen, og manglende ferdigstilte resultater fra
klassifikatoren. 

Kanskje det første neste trinnet er å presisere det vitenskapelige
målet med ditto  evaluering. Det er ikke noe
åpenbart behov for et mer komplisert nevralt nettverk eller en
omfattende ny serie DFT-beregninger før dette er gjort. En grundig
undersøkelse som viser begrensningene ved overføring fra et lite og
skjevt sammensatt datasett av vertsmaterialer, kan i seg selv være
vitenskapelig nyttig, forutsatt at analysen er fullført og
konklusjonene følger av resultatene.

Noe jeg syns er bra er det du har om dkillet mellom defektbasert
emisjon og emisjon fra kvanteprikker eller andre nanostrukturer,
fysisk fornuftig. Skillet mellom utforskende forbehandling av hele
datasettet og forbehandling som bare tilpasses treningsdataene, er
også et godt valg: På s. 33 står det uttrykkelig at klassifikatorens
verktøy for imputering, variansfiltrering og skalering tilpasses
treningsdelsettet. Dette bør videreføres og kontrolleres i
implementasjonen. PCA er et rimelig verktøy for utforskende analyse,
og diskusjonen av manglende data viser oppmerksomhet om et viktig
praktisk problem. Oppgaven bør bygge videre på disse styrkene kanskje
avgrense andre.

Men en avklaring, vva skal predikeres? Det du sendte veksler mellom at
et materiale er en enkeltfotonemitter, at det kan være vert for en
egnet defekt, at det egner seg for en spinnkvantebit, og at det ligner
kjente vertsmaterialer. Dette er forskjellige mål. Et
bulkvertsmateriale er ikke i seg selv en spesifisert emitter. Det
fysiske systemet omfatter krystallfasen, defekten eller urenheten,
ladningstilstanden, det lokale miljøet og betingelsene for
fremstilling og måling.

Et  hovedspørsmål kan være: Kan deskriptorer for
sammensetning, struktur og elektroniske egenskaper i bulkmaterialer
brukes til å prioritere krystallfaser som ligner etablerte
vertsmaterialer for defektbasert enkeltfotonemisjon, og hvor robust er
prioriteringen overfor valg av etiketter og tilbakeholding av hele
kjemiske familier? Dette samsvarer langt bedre med det utførte
arbeidet enn en generell påstand om å kunne predikere om et vilkårlig
materiale egner seg for enkeltfotonemisjon.

Innledningen er foreløpig tom, så kanskje denne denne avklaringen er noe du kan fokusere på?.
Dvs, beskriv målet, hvorfor screening av vertsmaterialer er nyttig,
hvilken informasjon inngangsdeskriptorene inneholder og mangler, de
konkrete forskningsspørsmålene og kandidatens selvstendige bidrag
sammenlignet med Hebnes og annet tidligere arbeid. Et enkelt diagram
over vertsmaterialets energibånd og defektnivåer vil bidra til å
knytte fysikken til maskinlæringsoppgaven.

Dokumentasjonen bak de positive etikettene bør gjennomgås før videre
tuning. Avsnitt 3.2.1 og 4.1.1 samt tabell 6.1 i vedlegget er
sentrale. En opptelling i vedlegget gir:

| Kjemisk formel | Positive databaseoppføringer |
|---|---:|
| SiC | 26 |
| C | 11 |
| Si | 4 |
| AlN | 2 |
| GaN | 2 |
| AlAs | 1 |
| Totalt | 46 |

De 46 positive oppføringene dekker altså bare seks kjemiske formler,
og 37 tilhører SiC eller C. De utgjør ikke 46 uavhengige
eksperimentelle demonstrasjoner fordelt på kjemisk forskjellige
vertsmaterialer. Teksten på s. 36 oppgir 12 C-oppføringer; de oppgitte
antallene summerer seg derfor til 47 i stedet for 46.

Mer grunnleggende er det slik at eksperimentell dokumentasjon av en
emitter i én fase ikke viser at alle polymorfer med samme kjemiske
formel egner seg som vertsmaterialer. At en oppføring i Materials
Project samsvarer med en eksperimentelt rapportert krystallstruktur,
dokumenterer strukturens opphav, ikke eksperimentell observasjon av
enkeltfotonemisjon i denne strukturen. Karbonallotroper og
SiC-polytyper gjør dette skillet særlig viktig. Utkastet dokumenterer
foreløpig ikke det fasespesifikke grunnlaget som trengs for å begrunne
den nåværende utvidelsen av etikettene.

For hver positiv oppføring bør det oppgis database-ID, kjemisk formel,
romgruppe eller polytype, eksperimentell fase, emisjonsmekanisme,
defektens identitet når den er kjent, relevant eksperimentell
dokumentasjon, primærreferanse og en kategori for hvor sikker
etiketten er. Skill mellom en fase som er direkte dokumentert, og en
beslektet fase som har fått en foreløpig etikett. Det er akseptabelt
at datasettet med sikre etiketter blir mindre etter denne
gjennomgangen.

Behold fysisk forskjellige polymorfer som egne strukturer, men ikke
anta at de har samme optiske oppførsel. Bruk struktursammenligning til
å identifisere duplikater eller nesten identiske oppføringer, og hold
beslektede oppføringer samlet under valideringen. En følsomhetsanalyse
med én representant per kjemisk formel, eller jevnere vekting av
familiene, vil vise hvor sterkt overrepresentasjonen av SiC og karbon
påvirker analysen.


Kan du sjekke om 
begrunnelsen hentet fra Ferrenti gjelder ønskelige egenskaper for
bestemte anvendelser av atomlignende kvantedefekter, blant annet
spinnkoherens og optisk koherens. Kan du sjekke? Klarte det ikke sjøl.


Har du mulighet til å undersøke om magnetisme, polaritet,
symmetrietiketter eller nært beslektede indikatorer inngår blant
modellens kjennetegn? Hvis en negativ klasse defineres ved en egenskap
som modellen mottar direkte, kan høy klassifikasjonsnøyaktighet ganske
enkelt gjenspeile at modellen gjenfinner utvalgsregelen. Sammenlign
resultatene med og uten disse kjennetegnene før det hevdes at modellen
tilfører prediktiv informasjon. Manuskriptet gir ikke nok detaljer til
å fastslå at denne bestemte snarveien forekommer i implementasjonen.

En tilnærming med positive og umerkede eksempler, eller en
enklassetilnærming, er verdt å vurdere, men er ingen automatisk
løsning. Blant de kjente positive eksemplene er godt studerte kjemiske
systemer sterkt overrepresentert, og mange metoder for læring fra
positive og umerkede eksempler bygger på antakelser om hvordan de
positive er valgt ut. For en håndterbar oppgave kan en uttrykkelig
avgrenset studie av likhetsbasert rangering være tilstrekkelig; det er
valgfritt å implementere flere læringsparadigmer.

**Fullfør enkle referansemodeller før søket etter nevrale nettverk
  utvides.** En liten sammenligning er tilstrekkelig: en konstant
  eller enkel standardprediktor, en transparent referanse basert på
  sammensetning eller likhet, en regularisert logistisk klassifikator
  og ett kompakt nevralt nettverk. Bruk de samme ytre oppdelingene,
  den samme informasjonen og de samme evalueringsmålene for alle. En
  ekstra trebasert modell er valgfri. En sterk enkel referansemodell
  vil være et nyttig resultat, ikke et problem for oppgaven.

Søket på s. 34 med 100 hyperparameterforsøk på et svært lite
valideringssett kan selektere tilfeldige mønstre i dette settet. Søket
oppgir dessuten dybde, men spesifiserer ikke fullt ut bredden på de
skjulte lagene, batchstørrelse, initialisering, endelig arkitektur
eller endelig treningsprosedyre. Hold søket begrenset, dokumenter det,
og foreta modellvalget innenfor treningsdataene, med gruppert indre
validering der det er gjennomførbart. Hvis dette blir for ustabilt med
de tilgjengelige gruppene, velg noen få enkle konfigurasjoner på
forhånd og oppgi begrensningen. En ytre evaluering bør ikke brukes
gjentatte ganger til å endre modellen.

Rapporter den faktiske endelige arkitekturen og antallet parametere,
målfunksjon og evalueringsmål, trenings- og valideringskurver,
stoppkriterium og resultater på tvers av relevante oppdelinger eller
tilfeldige startverdier. Balansert nøyaktighet, klassevis
gjenfinningsgrad og presisjon, en forvekslingsmatrise og et
rangeringsmål er rimelige valg, avhengig av oppgaven. Vis antall i
tillegg til prosentandeler; i et testsett med omtrent ti oppføringer
endres nøyaktigheten med rundt ti prosentpoeng for hver prediksjon som
endrer korrekthet. Gjentatte kjøringer med ulike startverdier skaper
ikke nye uavhengige kjemiske eksempler.

Ikke tolk en sigmoidutgang som sannsynligheten for at en kandidat
eksperimentelt vil kunne være vert for en
enkeltfotonemitter. Treningen bruker kunstig utvalgte, nesten
balanserte klasser og usikre negative eksempler; forekomsten og
observasjonsprosessen i hele databasen er ukjent. En modellskår kan
brukes til rangering uten å ha denne fysiske sannsynlighetstolkningen.


Begrunnelsen for å velge min-maks-skalering, nemlig at de fleste
kjennetegnene ikke er normalfordelte, er
utilstrekkelig. Standardisering krever ikke normalfordeling. Velg og
vurder skalering ut fra deskriptorenes enheter, uteliggere og
modellens geometri. Medianimputering og min-maks-skalering er mulige
valg, men virkningene bør undersøkes fremfor å fremstilles som entydig
påkrevd av fordelingene.

Det tydelige visuelle skillet mellom 2D- og 3D-dataene kan skyldes
reell fysikk knyttet til dimensjonalitet, effekter av databasekilde,
manglende data eller artefakter i representasjonen. De nåværende
figurene kan ikke skille mellom disse forklaringene. For eksempel
bruker 2DMatPedia en periodisk supercelle med et betydelig
vakuumområde. Volumavhengige deskriptorer må derfor behandles med
særlig omtanke ved sammenligning med bulkdeskriptorer. Dette er en
bekymring som kan testes, ikke en dokumentert forklaring på
figurene. Den opprinnelige databaseartikkelen angir minst 20 Å
vakuumseparasjon,se https://www.nature.com/articles/s41597-019-0097-3.

Slik treningen er beskrevet, har klassifikatoren ingen 2D-eksempler i
treningssettet. Skårene for 2D-materialer bør derfor behandles som
ekstrapolasjon. Den enkleste faglig forsvarlige avgrensningen er å
fullføre analysen av bulkvertsmaterialer og beholde 2D-materialene som
en separat utforskende sammenligning. En felles prediksjonsmodell vil
kreve en sterkere begrunnelse for at deskriptorer, etiketter og
validering er sammenlignbare.

**PCA er nyttig her, men den matematiske forklaringen må rettes.** På
  s. 22-23 er dimensjonene for observasjoner og kjennetegn byttet om,
  normaliseringen av kovariansen mangler, og det brukes en umulig
  identitet for en redusert basis. En konsistent konvensjon er en
  sentrert matrise X med dimensjon n_samples ganger p_features, og
  beholdte akser W med dimensjon p_features ganger k. Da er

\[
W^TW=I_k,\qquad Z=XW,\qquad \widehat X=ZW^T=XWW^T.
\]

Når k er mindre enn p_features, er WW^T en projektor, ikke
identitetsmatrisen. Rekonstruksjonen er generelt tilnærmet. For
usentrerte observasjoner må det tilpassede gjennomsnittet trekkes fra
før projeksjonen og legges til igjen etterpå. Utvalgskovariansen er
X^TX/(n_samples-1). PCA-komponentene er egenvektorer til denne
kovariansmatrisen, eller ekvivalent høyre singulærvektorer til
sentrert X. Det vanlige kvadratiske rekonstruksjonsmålet bruker
Frobenius-normen. Manuskriptets uttrykk som inneholder D D^T C, er
også dimensjonsmessig inkonsistent med den angitte matrisen av
komponentskårer.

Diskusjonen av sentroider på s. 29-31 blander sammen
enhetsnormaliserte akser og variansen til skårene langs disse
aksene. Vanlig PCA krymper ikke egenvektorene kunstig fordi de
forklarer mindre varians. Skårene har varianser gitt av de tilsvarende
egenverdiene. Hvis man ønsker en standardisert forskyvning mellom
gjennomsnitt, må den defineres eksplisitt, for eksempel

\[
d_j=\frac{\bar z_{A,j}-\bar z_{B,j}}{\sqrt{\lambda_j}}.
\]

Dette er en beskrivende forskyvning relativt til en valgt
referansevarians, ikke automatisk en signifikanstest. Å dele på
variansen i stedet for kvadratroten av den definerer en annen
størrelse. Retninger med små egenverdier kan forsterke støy, så disse
bør utelates eller regulariseres, og stabiliteten bør undersøkes. De
relevante skillene i programvaren mellom akser, forklart varians og
valgfri «whitening» er dokumentert i [PCA-referansen til
scikit-learn](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html).

**Eksperimentet med rekonstruksjonsfeil validerer foreløpig ikke
  prediksjon.** Avsnitt 4.1.2.2, s. 40, oppgir at alle de 46 positive
  observasjonene som PCA ble tilpasset til, har feil under 2 x
  10^-30. Et sentrert datasett med 46 observasjoner har rang på høyst
  45. Hvis alle komponentene som ikke er null, beholdes, vil
  datasettets egne rader per konstruksjon rekonstrueres med numerisk
  presisjon. De rapporterte feilene nær null er derfor ikke belegg for
  at metoden predikerer egnethet som vertsmateriale for en
  emitter. Oppgi antallet beholdte komponenter og den numeriske
  rangen, slik at beregningen kan tolkes.

Den foreslåtte skåren er likevel nyttig som et mål på avstand fra et
lært affint underrom. Begrensningen er lett å se: Et punkt på formen
mu + t w_1 har null residual for enhver t hvis w_1 er en beholdt akse,
selv om punktet ligger vilkårlig langt fra klyngen av treningsdata. En
lav residual dokumenterer verken lokal likhet, høy
sannsynlighetstetthet eller gunstig defektfysikk.

Bruk tilbakeholdte positive observasjoner til å evaluere skåren, helst
ved å holde tilbake hele familier med samme kjemiske formel. Rapporter
både residualavstanden fra underrommet og avstanden innenfor det
beholdte skårrommet, eller en enkel avstand til nærmeste nabo og en
markering av modellens gyldighetsområde. Velg antallet komponenter ut
fra treningsdataene, og undersøk følsomheten for rang,
deskriptorblokker og forbehandling. Et tilfeldig referanseunderrom
eller en enkel rangering etter avstand i sammensetning vil bidra til å
avklare om PCA tilfører nyttig informasjon. Presenter resultatet som
en utforskende likhets- eller avviksskår med mindre sterkere
validering gir grunnlag for mer.

Den rapporterte listen over de 70 øverst rangerte oppføringene er
opplysende. Den inneholder de 46 eksemplene som ble brukt i
tilpasningen. De resterende oppføringene beskrives som ni C-, sju Si-,
én SiC-, én AlAs-, én AlN- og én GaN-oppføring, samt fire ternære
ID-er. Dermed tilhører 20 av de resterende 24 oppføringene formler som
allerede er representert. De fire ternære ID-ene representerer bare to
ytterligere formler: AlGa3N4 og Ga3BN4. Disse er kandidater sett i
forhold til dette treningssettet, ikke etablerte nye emittere eller
dokumenterte oppdagelser. Det gjør ikke rangeringen ubrukelig; det
avgrenser hva resultatet foreløpig viser.

**Den øvrige statistikken bør støtte argumentasjonen, ikke bare gi
  flere beskrivende figurer.** Flere konkrete rettelser er nødvendige:

| Sted | Problem | Passende rettelse |
|---|---|---|
| s. 17, ligning (2.24) | Den oppgitte dimensjonen til vektmatrisen passer ikke med multiplikasjon av en kolonnevektor i R^n for å gi en vektor i R^m. | For h_next = a(W h + b) må W ha dimensjon m ganger n og b ha lengde m. |
| s. 19, ligning (2.27) | Binær kryssentropi blander sammen symbolene for prediksjon og målverdi og bruker et udefinert symbol. | Skriv -mean[y log(p) + (1-y) log(1-p)], med målverdi y i {0,1} og predikert p i (0,1). Kontroller implementeringen separat. |
| s. 25, ligning (2.36) | Cosinuslikhet oppgis å ligge i [0,1]. | For generelle reelle vektorer som ikke er null, er verdiområdet [-1,1]; PCA-akser kan ha begge fortegn. |
| s. 25, 37-39 | Rå cosinuslikhet mellom individuelle prinsipalkomponenter kan bli behandlet som strukturell ekvivalens. | Fortegnene til prinsipalkomponentene er vilkårlige; nesten degenererte akser kan rotere. Sammenlign absoluttverdien av orienteringstilpassede overlapp eller, helst, prinsipale vinkler mellom underrom, med samsvarende kjennetegnskoordinater og forbehandling. |
| s. 23-24 | LDA beskrives hovedsakelig som maksimering av avstanden mellom klassegjennomsnitt. | Bruk forholdet mellom spredning mellom klassene og spredning innenfor klassene; klasseseparasjon i dataene metoden er tilpasset til, validerer ikke de fysiske etikettene. |
| s. 23-24 og metodekapitlet | Langt flere kjennetegn enn merkede observasjoner gjør vanlig invertering av kovariansmatrisen problematisk. | Oppgi LDA-løseren og eventuell krymping av kovariansestimatet eller dimensjonsreduksjon basert bare på treningsdata; undersøk oppførselen på tilbakeholdte data. |
| s. 24-25, ligning (2.35) | ANOVA bruker et uferdig variansforhold og hevder at F > 1 betyr en meningsfull forskjell. | Bruk kvadratsummer dividert på tilhørende frihetsgrader og riktig referansefordeling; rapporter effektstørrelser og usikkerhet, og håndter multippel testing. F > 1 er ikke et signifikanskriterium. |

ANOVAs betingelse om lik varians gjelder gruppene for
responsvariabelen som testes, ikke at alle kjennetegn skal ha samme
varians. Avhengighet mellom polymorfer har også betydning. Vanlige
parametriske p-verdier bør ikke tolkes som om hver beslektede struktur
var et uavhengig eksperiment.

PCA fanger opp varians, ikke nødvendigvis informasjon som er nyttig
for prediksjon. LDA er uttrykkelig konstruert for å finne separasjon
mellom etikettene. Ingen av metodene beviser at etikettene har en
uavhengig korrekt fysisk betydning. Før enda en statistisk metode
legges til, bør det avklares hvilket konkret spørsmål den
besvarer. Det vil være rimelig å beholde PCA, en tydelig likhetsskår
og en liten sammenligning av prediktive modeller, og flytte sekundære
analyser til et vedlegg eller utelate dem.

Den skrevne DFT-bakgrunnen inneholder vesentlige, men korrigerbare
  feil.** Disse bør rettes før innlevering, samtidig som kapitlet
  konsentreres om databaseegenskapene som brukes i prosjektet.

På s. 11, ligning (2.3), har den fulle Hamilton-operatoren kinetiske
ledd med positive Laplace-operatorer, ekstra faktorer på en halv og en
feilaktig halvert tiltrekning mellom elektroner og kjerner. Summer
over partikler av samme type må også utelukke par der en partikkel
telles sammen med seg selv, og ulike kjernetyper må ha forskjellige
ladninger. Et konsistent SI-uttrykk er

\[
\hat H=-\sum_A\frac{\hbar^2}{2M_A}\nabla_A^2
-\sum_i\frac{\hbar^2}{2m_e}\nabla_i^2
+\frac{e^2}{4\pi\epsilon_0}\left[
\sum_{A<B}\frac{Z_AZ_B}{|\mathbf R_A-\mathbf R_B|}
+\sum_{i<j}\frac{1}{|\mathbf r_i-\mathbf r_j|}
-\sum_{i,A}\frac{Z_A}{|\mathbf r_i-\mathbf R_A|}\right].
\]

Atomære enheter vil være like akseptabelt dersom dette oppgis og
brukes konsekvent. Ligning (2.5) blander nå dimensjonerte kinetiske
prefaktorer med et Coulomb-ledd i atomære enheter.

På s. 14 legges V_aux til etter at T_aux er trukket fra i definisjonen
av utvekslings- og korrelasjonsenergi rett før ligning (2.19). Erstatt
den tvetydige notasjonen med standardoppspaltingen, i atomære enheter,

\[
E[n]=T_s[n]+\int n(\mathbf r)v_{\rm ext}(\mathbf r)d^3r
+E_H[n]+E_{xc}[n],
\]

\[
E_H[n]=\frac12\iint\frac{n(\mathbf r)n(\mathbf r')}{|\mathbf r-\mathbf r'|}d^3r\,d^3r',
\qquad E_{xc}[n]=(T[n]-T_s[n])+(V_{ee}[n]-E_H[n]).
\]

Forklar at det ikke-vekselvirkende hjelpesystemet gjengir
tettheten. Dette innebærer ikke at den sanne vekselvirkende
bølgefunksjonen er et produkt av uavhengige
ett-elektronfunksjoner. Definer V i ligning (2.20) eksplisitt som
summen av det ytre potensialet og Hartree-potensialet dersom det er
dette som menes. Det sammenslåtte symbolet alene dokumenterer ikke at
Hartree-leddet er utelatt.

På s. 15 skriver du   at selvkonsistens mellom påfølgende tettheter
innebærer at grunntilstanden er funnet. Kan du sjekke det, jeg er ikke sikker. Det gir en
selvkonsistent løsning for det valgte funksjonalet og det numeriske
oppsettet. Konkurrerende okkupasjoner eller magnetiske konfigurasjoner
kan kreve energisammenligninger. Numeriske toleranser og
tetthetsblanding bør nevnes kort. På s. 13 er heller ikke
stasjonaritet under en normaliseringsbetingelse tilstrekkelig til å
identifisere grunntilstanden; eksiterte egentilstander er også
stasjonære. Bruk normaliserte tillatte bølgefunksjoner eller
Rayleigh-kvotienten, med riktig integrasjonsmål for mange elektroner.

På s. 12 og 15 krever orbitalenes ortonormalitet et indreprodukt, ikke
den punktvise ligningen psi_i* psi_j = delta_ij. Skriv integralet over
rom og spinn. På s. 12 bør utveksling og korrelasjon ikke omtales som
ytterligere ikke-Coulombske vekselvirkninger. Hartree-Fock tar hensyn
til antisymmetri og utveksling innenfor en determinant og utelater
korrelasjon utover denne tilnærmingen.

Den grunnleggende Hohenberg-Kohn-diskusjonen er et rimelig
utgangspunkt. Presiser grunntilstandstetthet, fast partikkeltall og
tillatte tettheter. Tettheten integrert over rommet gir
elektronantallet; den er ikke en sannsynlighetsfordeling for ett
elektron normalisert til én. Okkupasjonstall bør defineres dersom
romorbitaler brukes i stedet for spinnorbitaler.

Det viktigste som mangler om DFT, er det uferdige avsnittet om
begrensninger på s. 15. Forklar hvordan tilnærmede funksjonaler for
utvekslings- og korrelasjonsenergi, numeriske valg og dataenes opphav
påvirker deskriptorene. En vanlig feil i båndgap fra semilokale
funksjonaler gjør ikke den valgte søketerskelen på 0,1 eV til et
fysisk kriterium for enkeltfotonemisjon. DFT-deskriptorer for
bulkmaterialet i grunntilstanden og en skalar grunntilstandsenergi gir
ikke direkte emisjonsspektre for defekter, stabilitet av
ladningstilstander, strålingseffektivitet eller
fotonstatistikk. Kohn-Sham-gap, fundamentale gap og optiske gap må
skilles fra hverandre. Denne forbindelsen til den faktiske
prediksjonsoppgaven er viktigere enn en lang utledning av generell
kvantemekanikk. [Den opprinnelige screeningstudien til Weber
mfl.](https://qs3.mit.edu/images/pdf/Weber.pdf) gir nyttig kontekst
for skillet mellom vertsmaterialets egenskaper og defektspesifikke
krav.

**Innledningen om halvledere og emittere trenger også en kort
  gjennomgang av faglig presisjon.** På s. 9 er påstanden om at
  halvledere og isolatorer ikke kan lede strøm, for generell:
  Eksitasjon og doping skaper ladningsbærere, og hull bidrar også til
  strømledning. Metaller bør beskrives ut fra båndfylling eller
  båndoverlapp ved Fermi-nivået, fremfor en universell betegnelse på
  gapet. På s. 10 må en SPE defineres gjennom fotonstatistikk, ikke
  bare ved at den er en halvleder eller emitterer ved én
  bølgelengde. Oppgaven handler om faststoffemittere, som er en
  undergruppe av SPE-systemer. Lokaliserte tilstander i båndgapet
  dokumenterer ikke i seg selv en brukbar
  enkeltfotonovergang. Polarisasjon er én mulig koding av en
  fotonkvantebit, ikke en påstand om at fotoner bare har to optiske
  moder.

Behold bakgrunnen om spinnkvantebiter bare i den grad den trengs for
det faktiske målet og antakelsene bak screeningen. Hvis analysen
handler om vertsmaterialer for enkeltfotonemisjon og ikke om
spinnkvantebiters koherens, bør det sies tydelig. En kort og presis
fremstilling er å foretrekke fremfor tomme overskrifter som dekker et
bredere prosjekt enn det som er gjennomført.

Definer andreordens intensitetskorrelasjon g^(2)(tau), og forklar
g^(2)(0) < 1/2 som en mye brukt indikator på at emisjonen hovedsakelig
består av enkeltfotoner, med tilhørende måle- og
bakgrunnskontekst. Denne verdien alene beviser ikke entydig at det
finnes én mikroskopisk emitter. Skill mellom enkeltfotonrenhet,
lysstyrke, spektral koherens og spinnkoherens i stedet for å behandle
dem som om de var det samme. De trenger ikke alle å bli
prediksjonsmål; det er tilstrekkelig å identifisere hvilken egenskap
oppgaven undersøker.

**Hva resultatkapitlet foreløpig kan underbygge.** Det underbygger
  foreløpige observasjoner om geometrien til de valgte, forbehandlede
  deskriptorene og en liste over materialer som ligger nær
  PCA-underrommet tilpasset de positive observasjonene. Det
  underbygger foreløpig ikke en validert suksessrate for
  klassifikasjon, generalisering til ukjente kjemiske sammensetninger,
  kalibrerte sannsynligheter for egnethet eller nye emittere som
  fungerer eksperimentelt. Avsnitt 4.2.1 og 4.2.2 er i stor grad
  uskrevne, så modellens ytelse kan ikke vurderes ut fra denne PDF-en.

Flere numeriske og visuelle småting bør/må avklares med utgangspunkt i de faktiske resultatene fra den endelige analysen:

| Sted | Nødvendig rettelse |
|---|---|
| s. 35, tekst; s. 36, figur 4.1 | Figuren angir 2 prinsipalkomponenter = 27 %, 34 = 80 %, 75 = 90 % og 1 689 omtrent 100 %. Figurteksten oppgir andre antall og andeler; brødteksten sier også at de to første forklarer 10 %. Generer alle tre på nytt fra ett og samme resultat. |
| s. 30, figur 3.1 | Figurteksten beskriver et panel c som ikke finnes på den gjengitte siden. Rett figurteksten eller sett inn panelet igjen. |
| s. 29 og andre steder | Henvisninger til figur 4.5 stemmer ikke konsekvent med analysen som vises. Rett krysshenvisningene etter at resultatene er organisert. |
| s. 27-28, tabell 3.1 og databeskrivelsen | MP-antallet i tabellen ser ut til å bruke det samlede antallet, mens teksten oppgir 30 699 MP-oppføringer pluss 2 663 2D-oppføringer. Avklar totalantallet fra hver kilde og frafallet i hvert trinn. |
| s. 28, 33, 40 | Antall kjennetegn på 5 013, 4 975, 3 445 og 3 574 kan gjelde ulike trinn eller versjoner, men sammenhengen kan ikke reproduseres. Ta med en tabell som viser dataflyten. |
| s. 39, figur 4.4 | Etikettene står tett, og navn på kjennetegn er utelatt. Vis et mindre utvalg med fullstendige identifikatorer, og vis samlede bidrag fra deskriptorfamilier separat. Påstanden i teksten om at komponenten hovedsakelig består av Radialdistribution, må også kontrolleres mot den synlige deskriptorlisten. Tolk komponentladningene fremfor å behandle størrelsen deres som kausal betydning. |
| s. 41, figur 4.5 | PDF-en inneholder tomme LaTeX-bokser med filnavn i stedet for kandidatfigurene. De foreslåtte plasseringene av kandidatene kan derfor ikke kontrolleres visuelt. |
| s. 29, 31, 39-40 | Skill metodedefinisjoner fra empiriske resultater og tolkning. Fjern utkastkommentarer først etter at de underliggende spørsmålene er avklart. |

Unngå å bruke mye plass på å beskrive formen til hver synlige kurve
eller hvert enkeltstående punkt i to prinsipalkomponenter. Forklar hva
figuren undersøker, identifiser viktige avvikende punkter med
dokumentasjon om fase og deskriptorer, og oppgi begrensningene ved
projeksjonen. Tjuesju prosent forklart varians kan fortsatt være
visuelt nyttig, men kan ikke beskrive hele klassestrukturen. Omvendt
er en høy kumulativ forklart varians ikke belegg for at optiske
egenskaper kan predikeres.

**Gjør den korte kandidatlisten til en fysisk vurdering.** De fire
  rapporterte ternære ID-ene på s. 40 gir et håndterbart
  utgangspunkt. Vis for hver av dem den nøyaktige fasen og strukturen,
  dataenes opphav, bulkmaterialets båndgap og beregningsmetoden,
  relativ stabilitet, de nærmeste kjente vertsmaterialene og hvor
  stabil rangeringen er på tvers av rimelige analysevalg. Forklar om
  likheten ganske enkelt gjenspeiler en legering av
  nitridsammensetninger som allerede er representert. Undersøk
  eksisterende litteratur før kandidaten omtales som ny.

Identifiser hva som fortsatt må til for å dokumentere en emitter: en
plausibel lokalisert defekt eller urenhet, mulige ladningstilstander,
egnede lokaliserte nivåer og optiske overganger, samt konkurrerende
ikke-strålende prosesser eller ioniseringsprosesser. Disse kravene
hører hjemme i tolkningen; de trenger ikke alle å bli nye beregninger
i denne masteroppgaven. En grundig litteratur- og databasevurdering av
noen få kandidater er et realistisk neste steg. En liten, validert
defektberegning kan styrke oppgava men den bør ikke fortrenge
retting av etiketter og grundig validering.

**En håndterbar arbeidsrekkefølge.** Følgende rekkefølge er ment å bidra til å fullføre en sammenhengende oppgave uten å utvide prosjektet i det uendelige.

| Prioritet | Arbeid | Konkret kriterium for ferdigstillelse |
|---|---|---|
| 1 | Avklar mål og etikettenes dokumentasjonsgrunnlag | Ett eksplisitt hovedspørsmål; en dokumentasjonstabell for positive eksempler; terminologi som synliggjør at de negative etikettene er indirekte; kontrollerte overlapp mellom etiketter og kontrollerte fasetilordninger. |
| 2 | Etabler evalueringsopplegget | Materiallister der duplikater er fjernet og beslektede materialer er gruppert, heltallsantall i hvert delsett, forbehandling basert bare på treningssettet og et definert eksperiment der hele familier holdes tilbake. |
| 3 | Fullfør en liten modellsammenligning | Enkle referansemodeller og ett kompakt nevralt nettverk evaluert på de samme delsettene, med antall, kurver, usikkerhet/variasjon og dokumentert modellvalg. |
| 4 | Gjør PCA-påstandene faglig holdbare | Korrekte formler, oppgitt numerisk rang, gjenfinning av tilbakeholdte eksempler, residual og avstand innenfor underrommet, samt en begrenset følsomhetsstudie av forbehandling og deskriptorer. |
| 5 | Tolk noen få kandidater fysisk | En kort tabell og nærmere gjennomgang av enkeltkandidater som skiller mellom likhet, plausibilitet som vertsmateriale og manglende dokumentasjon om defekter. |
| 6 | Fullfør den vitenskapelige fremstillingen | Presis og konsentrert bakgrunn; ferdige resultater og diskusjon; innledning, sammendrag og konklusjon som samsvarer med det dokumenterte omfanget; rettede figurer og referanser. |

