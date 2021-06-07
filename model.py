import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

genotypes = ['RCH','RCh','RcH','Rch',
             'rCH','rCh','rcH','rch']
environments = ['neutral','bacteriophage','antibiotic']

def makefitmat(sp,sm):
    # make the fitness matrix given parameters
    # fm[environment][genotype]
    return [[1-2*sm,1-sm,1-sm,1,
             1-2*sm,1-sm,1-sm,1],
            [(1+sp)*(1-2*sm),(1+sp)*(1-sm),1-sm,1,
             (1+sp)*(1-2*sm),(1+sp)*(1-sm),1-sm,1],
            [(1+sp)*(1-2*sm),(1+sp)*(1-sm),(1+sp)*(1-sm),1+sp,
             1-2*sm,1-sm,1-sm,1]
           ]

def cg(genotype,allele):
    # check if allele is in the genotype
    if allele in genotype:
        return True
    else:
        return False

def h(g1,g2,gh,gH):
    # calcuate the probability of gene transfer success given alleles
    if cg(g1,'H'):
        if cg(g2,'H'):
            return 2*gH
        else:
            return gH
    else:
        if cg(g2,'H'):
            return gH
        else:
            return gh

def getEvents(envs):
    # get start/end of each threath even
    dif = [envs[i]-envs[i+1] for i in range(0,len(envs)-1)]
    starts = [x+1 for x,d in enumerate(dif) if d < 0]
    ends = [x+1 for x,d in enumerate(dif) if d > 0]
    return starts,ends


def env_model(env,n,l,t):
    # model environmental states in the model
    envs = [0]*n
    event = np.random.randint(0,(n-l)/2)
    if env == 'Singular':
        envs[event:event+l] = [t]*l
    elif env == 'Cyclical':
        while event+2*l < n:
            envs[event:event+l] = [t]*l
            event += 2*l
    elif env == 'Alternating':
        while event+2*l < n:
            envs[event:event+l] = [t]*l
            event += 2*l
            if t == 2:
                t = 1
            else:
                t = 2
    elif env == 'Neutral':
        envs = [0]*n
    return envs

def transfer(freqs,gh,gH):
    # transfer step frequencies
    transfered = []
    for g,f in enumerate(freqs):
        successes = [h(genotypes[g],genotypes[i],gh,gH) for i in range(0,4)]
        transfers = sum([f*freqs[i]*s for i,s in enumerate(successes)])
        # print(genotypes[g],successes,transfers)
        if cg(genotypes[g],'R'):
            ft = f + transfers
        else:
            ft = f - transfers
        transfered.append(ft)
    return transfered

def mutation(freqs,env,mu_r,mu_R):
    # mutation step frequencies
    mutated = []
    for g,f in enumerate(freqs):
        ng = freqs[(g+4)%len(genotypes)]
        if cg(genotypes[g],'R'):
            if env == 1: # antibiotic env
                mu = np.sqrt(mu_R)
            else:
                mu = mu_R
        else:
            if env == 1: # antibiotic env
                mu = np.sqrt(mu_r)
            else:
                mu = mu_r
        tt = (1-mu)*f+mu*ng
        mutated.append(tt)
    return mutated

def selection(freqs,env,fitmat):
    # selection step frequencies
    selected = []
    avg_fitness = sum([f*fitmat[env][i] for i,f in enumerate(freqs)])
    for g,f in enumerate(freqs):
        st = f*fitmat[env][g]/avg_fitness
        selected.append(st)
    return selected

def simulate(init,n,envs,gh,gH,mu_r,mu_R,fitmat):
    # run the simulation for n generations given params
    generations = [init]
    for i in range(n):
        last_gen = generations[-1]
        transfered = transfer(last_gen,gh,gH)
        mutated = mutation(transfered,envs[i],mu_r,mu_R)
        selected = selection(mutated,envs[i],fitmat)
        generations.append(selected)
    return generations

def plotSim(generations,envs,title):
    # plot a single simulation run
    starts, ends = getEvents(envs)
    for s,e in zip(starts,ends):
        plt.axvline(x=s,color='g',linestyle='--')
        plt.axvline(x=e,color='r',linestyle='--')

    x = list(range(0,n+1))
    geno_frequencies = list(map(list, zip(*generations)))
    for genotype,freqs in enumerate(geno_frequencies):
        plt.plot(x,freqs,label=genotypes[genotype])
    plt.title(title)
    plt.ylabel('Frequency')
    plt.xlabel('Generation')
    plt.legend()
    plt.tight_layout()
    path = './figures/{}'.format(title.replace(' ','_'))
    plt.savefig(path)
    plt.close()

def plotScenario(init,gh,gH,mu_r,mu_R,fitmat,n,scenario,l):
    # plot the antibitoic and phage threats for the specificed environmental scenario
    envs = env_model(scenario,n,l,1)
    generations = simulate(init,n,envs,gh,gH,mu_r,mu_R,fitmat)
    plotSim(generations,envs,'{} Phage'.format(scenario))
    if scenario ==  'Alternating':
        envs = [x if x != 1 else -1 for x in envs]
        envs = [x if x != 2 else 1  for x in envs]
        envs = [x if x != -1 else 2 for x in envs]
    else:
        envs = [x if x != 1 else 2 for x in envs]
    generations = simulate(init,n,envs,gh,gH,mu_r,mu_R,fitmat)
    plotSim(generations,envs,'{} Antibiotic'.format(scenario))

def plotTerminalFrequencyStability(init,gh,gH,mu_r,mu_R,fitmat,n):
    # plot final generation frequencies over various environmental turnover rates (l)
    stabilities = range(10,100)
    final_neutral, final_antibiotic, final_phage, final_alternating = [], [], [], []

    for l in tqdm(stabilities,desc='stability...'):
        # No Events
        envs = env_model('Neutral',n,l,0)
        generations = simulate(init,n,envs,gh,gH,mu_r,mu_R,fitmat)
        final_neutral.append(generations[-1])
        # Cyclical Phage
        envs = env_model('Cyclical',n,l,1)
        generations = simulate(init,n,envs,gh,gH,mu_r,mu_R,fitmat)
        final_phage.append(generations[-1])
        # Cyclical Antibiotic
        envs = env_model('Cyclical',n,l,2)
        generations = simulate(init,n,envs,gh,gH,mu_r,mu_R,fitmat)
        final_antibiotic.append(generations[-1])
        # Alternating
        envs = env_model('Alternating',n,l,1)
        generations = simulate(init,n,envs,gh,gH,mu_r,mu_R,fitmat)
        final_alternating.append(generations[-1])

    fig, axs = plt.subplots(2,2,sharex=True,sharey=True,constrained_layout=True)
    axs[0,0].set_title('No Events')
    axs[0,0].set_ylabel('Frequency in Final Generation')
    for genotype,freqs in enumerate(list(map(list, zip(*final_neutral)))):
        axs[0,0].plot(stabilities,freqs,label=genotypes[genotype])
    axs[1,0].set_title('Cyclical Antibiotic')
    axs[1,0].set_ylabel('Frequency in Final Generation')
    axs[1,0].set_xlabel('Environmental Turnover')
    for genotype,freqs in enumerate(list(map(list, zip(*final_antibiotic)))):
        axs[1,0].plot(stabilities,freqs,label=genotypes[genotype])
    axs[0,1].set_title('Cyclical Phage')
    for genotype,freqs in enumerate(list(map(list, zip(*final_phage)))):
        axs[0,1].plot(stabilities,freqs,label=genotypes[genotype])
    axs[1,1].set_title('Alternating Events')
    axs[1,1].set_xlabel('Environmental Turnover')
    for genotype,freqs in enumerate(list(map(list, zip(*final_alternating)))):
        axs[1,1].plot(stabilities,freqs,label=genotypes[genotype])

    fig.suptitle('Genotype Frequencies vs Environmental Turover')
    axs[0,0].legend()
    plt.savefig('./figures/stability.png')
    plt.close()


if __name__ == '__main__':
    # initial genotype frequencies
    init = [0.001,0.025,0.025,0.01,
            0.01,0.01,0.01]
    init.append(1-sum(init))
    # transfer
    gh, gH = 0.000001, 0.00005
    # mutation
    mu_r, mu_R = 0.00001, 0.00001
    # fitness
    sp, sm = 0.2, 0.01
    fitmat = makefitmat(sp,sm)
    n = 1000
    # run singular simulation
    scenario, l = 'Singular', 50
    plotScenario(init,gh,gH,mu_r,mu_R,fitmat,n,scenario,l)
    # run cyclic simulation
    scenario, l = 'Cyclical',  100
    plotScenario(init,gh,gH,mu_r,mu_R,fitmat,n,scenario,l)
    # run alternating simulation
    scenario, l = 'Alternating',  100
    plotScenario(init,gh,gH,mu_r,mu_R,fitmat,n,scenario,l)
    # plot env stability
    plotTerminalFrequencyStability(init,gh,gH,mu_r,mu_R,fitmat,n)
