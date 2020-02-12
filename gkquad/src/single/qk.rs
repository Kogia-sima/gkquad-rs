#![allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::uninit_assumed_init
)]

use std::borrow::Borrow;
use std::ops::{Add, AddAssign};

use super::common::{Integrand, Interval};
use super::qk_impl::qk;
use super::util::Aligned;

/// holds the result of Gauss-Kronrod integration
#[derive(Debug)]
pub struct QKResult {
    /// approximation to the integral
    pub estimate: f64,
    /// estimate of the modulus of the absolute error
    pub delta: f64,
    /// approximation to the integral of |f|
    pub absvalue: f64,
    /// approximation to the integral of |f - quad(f)/(b - a)|
    pub asc: f64,
}

impl<T: Borrow<QKResult>> Add<T> for QKResult {
    type Output = QKResult;

    #[inline]
    fn add(self, other: T) -> QKResult {
        let other = other.borrow();
        QKResult {
            estimate: self.estimate + other.estimate,
            delta: self.delta + other.delta,
            absvalue: self.absvalue + other.absvalue,
            asc: self.asc + other.asc,
        }
    }
}

impl<T: Borrow<QKResult>> AddAssign<T> for QKResult {
    #[inline]
    fn add_assign(&mut self, other: T) {
        let other = other.borrow();
        self.estimate += other.estimate;
        self.delta += other.delta;
        self.absvalue += other.absvalue;
        self.asc += other.asc;
    }
}

/// Performs Gauss-Kronrod integration with 17-point kronrod rule
#[inline]
pub fn qk17<F: Integrand>(f: &mut F, r: &Interval) -> QKResult {
    unsafe {
        let mut fv = Aligned::<[f64; 17]>::uninit();
        qk(f, r, &XGK17, &WG17, &WGK17, WCK17, &mut *fv)
    }
}

/// Performs Gauss-Kronrod integration with 25-point kronrod rule
#[inline]
pub fn qk25<F: Integrand>(f: &mut F, r: &Interval) -> QKResult {
    unsafe {
        let mut fv = Aligned::<[f64; 25]>::uninit();
        qk(f, r, &XGK25, &WG25, &WGK25, WCK25, &mut *fv)
    }
}

/// Performs Gauss-Kronrod integration with 33-point kronrod rule
#[inline]
pub fn qk33<F: Integrand>(f: &mut F, r: &Interval) -> QKResult {
    unsafe {
        let mut fv = Aligned::<[f64; 33]>::uninit();
        qk(f, r, &XGK33, &WG33, &WGK33, WCK33, &mut *fv)
    }
}

/// Performs Gauss-Kronrod integration with 41-point kronrod rule
#[inline]
pub fn qk41<F: Integrand>(f: &mut F, r: &Interval) -> QKResult {
    unsafe {
        let mut fv = Aligned::<[f64; 41]>::uninit();
        qk(f, r, &XGK41, &WG41, &WGK41, WCK41, &mut *fv)
    }
}

/// Performs Gauss-Kronrod integration with 49-point kronrod rule
#[inline]
pub fn qk49<F: Integrand>(f: &mut F, r: &Interval) -> QKResult {
    unsafe {
        let mut fv = Aligned::<[f64; 49]>::uninit();
        qk(f, r, &XGK49, &WG49, &WGK49, WCK49, &mut *fv)
    }
}

/// Performs Gauss-Kronrod integration with 57-point kronrod rule
#[inline]
pub fn qk57<F: Integrand>(f: &mut F, r: &Interval) -> QKResult {
    unsafe {
        let mut fv = Aligned::<[f64; 57]>::uninit();
        qk(f, r, &XGK57, &WG57, &WGK57, WCK57, &mut *fv)
    }
}

// Gauss-Kronrod weights
// source: https://keisan.casio.com/exec/system/1289382036

const XGK17: [f64; 8] = [
    0.183434642495649804939476142360184,
    0.360701097928131957192548622296891,
    0.525532409916328985817739049189246,
    0.672354070945158677156310738092831,
    0.796666477413626739591553936475830,
    0.894120906847456421948361017538251,
    0.960289856497536231683560868569473,
    0.993379875881716155935888069019671,
];

const WG17: [f64; 4] = [
    0.362683783378361982965150449277196,
    0.313706645877887287337962201986601,
    0.222381034453374470544355994426241,
    0.101228536290376259152531354309962,
];

const WCK17: f64 = 0.184446405744691643528970955705643;

const WGK17: [f64; 8] = [
    0.181400025068034643061748525172550,
    0.172070608555211311857294880203857,
    0.156652606168188400490248088486969,
    0.136263109255172215262338745254506,
    0.111646370826839613222108158933941,
    0.082482298931358330688625193445608,
    0.049439395002139308500363969446998,
    0.017822383320710355152786961202750,
];

const XGK25: [f64; 12] = [
    0.125233408511468915472441369463853,
    0.248505748320469276267790960362718,
    0.367831498998180193752691536643718,
    0.481339450478157092935943615018832,
    0.587317954286617447296702418940534,
    0.684059895470055893944929100341155,
    0.769902674194304687036893833212818,
    0.843558124161153244792141885059839,
    0.904117256370474856678465866119096,
    0.950537795943121296549060195131619,
    0.981560634246719250690549090149281,
    0.996933922529595426912350237258385,
];

const WG25: [f64; 6] = [
    0.249147045813402785000562436042951,
    0.233492536538354808760849898924878,
    0.203167426723065921749064455809798,
    0.160078328543346226334652529543359,
    0.106939325995318430960254718193996,
    0.047175336386511827194615961485017,
];

const WCK25: f64 = 0.125556893905474335304296132860078;

const WGK25: [f64; 12] = [
    0.124584164536156073437312473209229,
    0.121626303523948383246099758091310,
    0.116712053501756826293580745305730,
    0.110022604977644072635907398742250,
    0.101649732279060277715688770491228,
    0.091549468295049210528171939739614,
    0.079920275333601701493392609529783,
    0.067250907050839930304940940047316,
    0.053697017607756251228889163320458,
    0.038915230469299477115089632285863,
    0.023036084038982232591084580367969,
    0.008257711433168395757693922439212,
];

const XGK33: [f64; 16] = [
    0.095012509837637440185319335424958,
    0.189168579018083726314712086634942,
    0.281603550779258913230460501460496,
    0.371483780878416288684027947079523,
    0.458016777657227386342419442983578,
    0.540407676352139743625428344051591,
    0.617876244402643748446671764048791,
    0.689741106681762303876199858016583,
    0.755404408355003033895101194847442,
    0.814240287062444468094572577588398,
    0.865631202387831743880467897712393,
    0.909157667012342948622557332882809,
    0.944575023073232576077988415534608,
    0.971505950969392594304109858259183,
    0.989400934991649932596154173450333,
    0.998239274145444514183282371262429,
];

const WG33: [f64; 8] = [
    0.189450610455068496285396723208283,
    0.182603415044923588866763667969220,
    0.169156519395002538189312079030360,
    0.149595988816576732081501730547479,
    0.124628971255533872052476282192016,
    0.095158511682492784809925107602246,
    0.062253523938647892862843836994378,
    0.027152459411754094851780572456018,
];

const WCK33: f64 = 0.095154216080498307020415055955840;

const WGK33: [f64; 16] = [
    0.094728401247230041326733968799471,
    0.093438674060921230478147190079644,
    0.091292032828191662272226594165447,
    0.088337502579112731413574030468030,
    0.084595803792590637585459382901613,
    0.080053941263719286445983306344537,
    0.074769823885599557533866452931194,
    0.068862995191531243001306696674629,
    0.062358806011834855470397014742265,
    0.055205633095422172303120001393629,
    0.047506215976407017234986327394423,
    0.039512951202421964003078606929551,
    0.031260543647380528239661154138738,
    0.022498859440049444029474671942737,
    0.013257930688091157245432985393750,
    0.004742777049247317906344087722423,
];

const XGK41: [f64; 20] = [
    0.076526521133497333754640409398838,
    0.152605465240922675505220241022678,
    0.227785851141645078080496195368575,
    0.301627868114913004320555356858592,
    0.373706088715419560672548177024927,
    0.443593175238725103199992213492640,
    0.510867001950827098004364050955251,
    0.575140446819710315342946036586425,
    0.636053680726515025452836696226286,
    0.693237656334751384805490711845932,
    0.746331906460150792614305070355642,
    0.795041428837551198350638833272788,
    0.839116971822218823394529061701521,
    0.878276811252281976077442995113078,
    0.912234428251325905867752441203298,
    0.940822633831754753519982722212443,
    0.963971927277913791267666131197277,
    0.981507877450250259193342994720217,
    0.993128599185094924786122388471320,
    0.998859031588277663838315576545863,
];

const WG41: [f64; 10] = [
    0.152753387130725850698084331955098,
    0.149172986472603746787828737001969,
    0.142096109318382051329298325067165,
    0.131688638449176626898494499748163,
    0.118194531961518417312377377711382,
    0.101930119817240435036750135480350,
    0.083276741576704748724758143222046,
    0.062672048334109063569506535187042,
    0.040601429800386941331039952274932,
    0.017614007139152118311861962351853,
];

const WCK41: f64 = 0.076600711917999656445049901530102;

const WGK41: [f64; 20] = [
    0.076377867672080736705502835038061,
    0.075704497684556674659542775376617,
    0.074582875400499188986581418362488,
    0.073030690332786667495189417658913,
    0.071054423553444068305790361723210,
    0.068648672928521619345623411885368,
    0.065834597133618422111563556969398,
    0.062653237554781168025870122174255,
    0.059111400880639572374967220648594,
    0.055195105348285994744832372419777,
    0.050944573923728691932707670050345,
    0.046434821867497674720231880926108,
    0.041668873327973686263788305936895,
    0.036600169758200798030557240707211,
    0.031287306777032798958543119323801,
    0.025882133604951158834505067096153,
    0.020388373461266523598010231432755,
    0.014626169256971252983787960308868,
    0.008600269855642942198661787950102,
    0.003073583718520531501218293246031,
];

const XGK49: [f64; 24] = [
    0.064056892862605626085043082624745,
    0.127851240286216699326667098746182,
    0.191118867473616309158639820757070,
    0.253600430369677910537489448545236,
    0.315042679696163374386793291319810,
    0.375191154797950870412699517733197,
    0.433793507626045138487084231913350,
    0.490612365463446416038345668979945,
    0.545421471388839535658375617218372,
    0.597990513906078379703258534153066,
    0.648093651936975569252495786910748,
    0.695532072396778802461733300689389,
    0.740124191578554364243828103099978,
    0.781677226476464317170710943246370,
    0.820001985973902921953949872669745,
    0.854953803904051366827147029446991,
    0.886415527004401034213154341982197,
    0.914240690794911500874994145088397,
    0.938274552002732758523649001708722,
    0.958441684405209407000907641216486,
    0.974728555971309498198391993008169,
    0.987040496015809045661879675246613,
    0.995187219997021360179997409700737,
    0.999201056021875051655759554001639,
];

const WG49: [f64; 12] = [
    0.127938195346752156974056165224695,
    0.125837456346828296121375382511184,
    0.121670472927803391204463153476262,
    0.115505668053725601353344483906784,
    0.107444270115965634782577342446606,
    0.097618652104113888269880664464247,
    0.086190161531953275917185202983743,
    0.073346481411080305734033615253117,
    0.059298584915436780746367758500109,
    0.044277438817419806168602748211338,
    0.028531388628933663181307815951878,
    0.012341229799987199546805667070037,
];

const WCK49: f64 = 0.064100463769266717756751728006586;

const WGK49: [f64; 24] = [
    0.063969626241376347509644034120127,
    0.063574878712972420734803194875307,
    0.062917112269698113170570246311836,
    0.062004001419781273939901272164948,
    0.060838034873127860983583676308816,
    0.059416543245952097694404582907668,
    0.057748675375690034539937419561017,
    0.055851715103063319555227688826167,
    0.053727945501751317395127210655595,
    0.051371951383457638487252496393224,
    0.048801386753259051450053450726567,
    0.046045891776563833411740269972464,
    0.043105928193695268092760628613405,
    0.039967623893622181222630944362560,
    0.036658100242212953886971295722458,
    0.033227378319829572434321287452346,
    0.029671306090659117491311160688321,
    0.025951194661374106405563815993065,
    0.022104084900061889640191567506721,
    0.018231178550387270198471476568716,
    0.014327444630883845848554484407082,
    0.010259786409280761245641108764037,
    0.006025671015720143722039610063186,
    0.002152308550946222060921213031076,
];

const XGK57: [f64; 28] = [
    0.055079289884034270426516527341880,
    0.109991649290732781027815261597935,
    0.164569282133380771281471777891166,
    0.218646765777540440144053724079084,
    0.272061627635178077676826356125770,
    0.324650637057820445835349056474133,
    0.376251516089078710221357209556087,
    0.426709291931107209807199466546646,
    0.475874224955118261034411847667434,
    0.523594465442833020107971645866323,
    0.569720471811401719308003283356431,
    0.614115617729884143662211574584255,
    0.656651094038864961219898176506743,
    0.697193171944897829505055730972558,
    0.735610878013631772028144510292534,
    0.771793447243433832029564566052428,
    0.805641370917179171447885955425278,
    0.837044208419339397285962669742181,
    0.865892522574395048942254567379687,
    0.892108688876811896164200241998283,
    0.915633026392132073869689423329927,
    0.936380950324850219172957334225778,
    0.954259280628938197254101839705216,
    0.969231347359594367303744365650804,
    0.981303165370872753694559945807830,
    0.990417243607036921336294784891355,
    0.996442497573954449950436390483311,
    0.999409527464458185430469787030878,
];

const WG57: [f64; 14] = [
    0.110047013016475196282376265601818,
    0.108711192258294135253571519303673,
    0.106055765922846417910416436996811,
    0.102112967578060769814216638505712,
    0.096930657997929915850489006095441,
    0.090571744393032840942186031336784,
    0.083113417228901218390396498244333,
    0.074646214234568779023931887173022,
    0.065272923966999595793397566775505,
    0.055107345675716745431482918226946,
    0.044272934759004227839587877653207,
    0.032901427782304379977630819170532,
    0.021132112592771259751500380993265,
    0.009124282593094517738816153922952,
];

const WCK57: f64 = 0.055107015299777746735626928438028;

const WGK57: [f64; 28] = [
    0.055023774809077271036502426453062,
    0.054772896844694716781166407068630,
    0.054354779876438122260022699069615,
    0.053773208209074282641162579514825,
    0.053029282478558707298691878688410,
    0.052121534240287052656082768398714,
    0.051054436965688590212331950140192,
    0.049836344907779101814673159246887,
    0.048468123831894560567113837884255,
    0.046947174743596009055304775696875,
    0.045282174971778382486767178658076,
    0.043486965463534754692944745996740,
    0.041561540422930381801750249005773,
    0.039500362038311021057620614398857,
    0.037316775066018058570339309681710,
    0.035032322351843966845355051650711,
    0.032644878891020445408809034380568,
    0.030141813915294599687594330310677,
    0.027542140838783441549627331827396,
    0.024881274360962674794539591845302,
    0.022153023381841145404571370933566,
    0.019325780376825872893481336671039,
    0.016424885799840213021325411478851,
    0.0135226653313631482127563325315156,
    0.010611858842171326615101945569492,
    0.007590801053104039604902428847390,
    0.004454936296692426768597724116182,
    0.001590736040706812893050065715678,
];
