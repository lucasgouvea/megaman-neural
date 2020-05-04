-- MegamanI/O by LucasGouvea
-- Intended for use with the BizHawk emulator and Megaman ROM.
-- mMke sure you have a save state named "DP1.state" at the beginning of a level,
-- and put a copy in both the Lua folder and the root directory of BizHawk.


buttonNames = {
	"B",
	"Left",
	"Down",
	"Right",
}

filename = "FREEZE.state"
scaleX = 890
scaleY = 2.76
frame = 0

-- Neural --

population = 50
deltaDisjoint = 2.0
deltaWeights = 0.4
deltaThreshold = 1.0
staleSpecies = 15
mutateConnectionsChance = 0.25
perturbChance = 0.90
crossoverChance = 0.75
linkMutationChance = 2.0
nodeMutationChance = 0.50
biasMutationChance = 0.40
stepSize = 0.1
disableMutationChance = 0.4
enableMutationChance = 0.2
timeoutConstant = 250
boxRadius = 6
maxNodes = 1000000
inputs = 2 + 1
outputs = #buttonNames


------------------------------------------------
------------------------------------------------
------------------------------------------------
----------------LIBS----------------------------
------------------------------------------------
------------------------------------------------

function clearJoypad()
	controller = {}
	for b = 1, #buttonNames do
		controller["P1 " .. buttonNames[b]] = false
	end
    joypad.set(controller)
end

function getPosition()
    megamanX = memory.read_u16_le(0XC04)
    megamanY = memory.readbyte(0xC08)
    megamanPos = {
        x = megamanX,
        y = megamanY
    }
    return megamanPos
end

function crossover(g1, g2)
	-- Make sure g1 is the higher fitness genome
	if g2.fitness > g1.fitness then
		tempg = g1
		g1 = g2
		g2 = tempg
	end

	local child = newGenome()
	
	local innovations2 = {}
	for i=1,#g2.genes do
		local gene = g2.genes[i]
		innovations2[gene.innovation] = gene
	end
	
	for i=1,#g1.genes do
		local gene1 = g1.genes[i]
		local gene2 = innovations2[gene1.innovation]
		if gene2 ~= nil and math.random(2) == 1 and gene2.enabled then
			table.insert(child.genes, copyGene(gene2))
		else
			table.insert(child.genes, copyGene(gene1))
		end
	end
	
	child.maxneuron = math.max(g1.maxneuron,g2.maxneuron)
	
	for mutation,rate in pairs(g1.mutationRates) do
		child.mutationRates[mutation] = rate
	end
	
	return child
end

function drawMap(genome)
    local network = genome.network
    local cells = {}
    local cell = {}
    local i = 1
    for dy = -boxRadius, boxRadius do
        for dx = -boxRadius, boxRadius do
            cell = {}
            cell.x = 50 + 5 * dx
            cell.y = 70 + 5 * dy
            if cell.y == 75 then
                cell.value = 1
            else
                cell.value = 0
            end                
            cells[i] = cell
            i = i + 1
        end
    end

    local biasCell = {}
    biasCell.x = 80
    biasCell.y = 110
    biasCell.value = 0
    cells[inputs] = biasCell


    gui.drawBox(49,71,51,78,0x00000000,0x80FF0000)
    for n, cell in pairs(cells) do
        if n > inputs or cell.value ~= 0 then
            local color = math.floor((cell.value + 1) / 2 * 256)
            if color > 255 then color = 255 end
            if color < 0 then color = 0 end
            local opacity = 0xFF000000
            if cell.value == 0 then
                opacity = 0x50000000
            end
            color = opacity + color * 0x10000 + color * 0x100 + color
            gui.drawBox(cell.x-2, cell.y-2, cell.x+2, cell.y+2, opacity, color)
        end
    end

    for o = 1, outputs do
		cell = {}
		cell.x = 220
		cell.y = 30 + 8 * o
		cell.value = network.neurons[maxNodes + o].value
		cells[maxNodes+o] = cell
		local color
		if cell.value > 0 then
			color = 0xFF0000FF
		else
			color = 0xFFFFFFFF
		end
		gui.drawText(188 , 24+14*o, buttonNames[o], color, 9)
	end

    
    color = 0xFF000000
end

function evaluateCurrent(megaX, megaY, bossX, bossY)
	local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]

    neuralInputs = getNeuralInputs(megaX, megaY, bossX, bossY)
	controller = evaluateNetwork(genome.network, neuralInputs)
	
	if controller["P1 Left"] and controller["P1 Right"] then
		controller["P1 Left"] = false
		controller["P1 Right"] = false
	end
	if controller["P1 Up"] and controller["P1 Down"] then
		controller["P1 Up"] = false
		controller["P1 Down"] = false
	end

	joypad.set(controller)
end

function getNeuralInputs(megaX, megaY, bossX, bossY)

    local neuralInputs = {}

    neuralInputs[1] = megaX - bossX
    neuralInputs[2] = megaY - bossY

	return neuralInputs
end

function initializePool(megaX, megaY, bossX, bossY)
	pool = newPool()

    for i = 1, population do
        basic = basicGenome()
		addToSpecies(basic)
	end

	initializeRun(megaX, megaY, bossX, bossY)
end

function initializeRun(megaX, megaY, bossX, bossY)
	savestate.load(filename);
	survivalMost = 0
	survival = 0
	pool.currentFrame = 0
	timeout = timeoutConstant
	clearJoypad()
	
	local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]
	generateNetwork(genome)
	evaluateCurrent(megaX, megaY, bossX, bossY)
end

function newNeuron()
	local neuron = {}
	neuron.incoming = {}
	neuron.value = 0.0
	
	return neuron
end

function newInnovation()
	pool.innovation = pool.innovation + 1
	return pool.innovation
end

function newPool()
	local pool = {}
	pool.species = {}
	pool.generation = 0
	pool.innovation = outputs
	pool.currentSpecies = 1
	pool.currentGenome = 1
	pool.currentFrame = 0
	pool.maxFitness = 0
	
	return pool
end

function containsLink(genes, link)
	for i = 1, #genes do
		local gene = genes[i]
		if gene.into == link.into and gene.out == link.out then
			return true
		end
	end
end

function newGene()
	local gene = {}
	gene.into = 0
	gene.out = 0
	gene.weight = 0.0
	gene.enabled = true
	gene.innovation = 0
	
	return gene
end

function pointMutate(genome)
	local step = genome.mutationRates["step"]
	for i = 1, #genome.genes do
		local gene = genome.genes[i]
		if math.random() < perturbChance then
			gene.weight = gene.weight + math.random() * step * 2 - step
		else
			gene.weight = math.random() * 4 - 2
		end
    end
    return genome
end

function linkMutate(genome, forceBias)
	local neuron1 = randomNeuron(genome.genes, false)
	local neuron2 = randomNeuron(genome.genes, true)
	 
	local newLink = newGene()
	if neuron1 <= inputs and neuron2 <= inputs then
		--Both input nodes
		return
    end

	if neuron2 <= inputs then
		-- Swap output and input
		local temp = neuron1
		neuron1 = neuron2
		neuron2 = temp
	end

	newLink.into = neuron1
	newLink.out = neuron2
	if forceBias then
		newLink.into = inputs
	end
	
	if containsLink(genome.genes, newLink) then
		return
    end

	newLink.innovation = newInnovation()
	newLink.weight = math.random()*4-2
	
	table.insert(genome.genes, newLink)
end

function randomNeuron(genes, nonInput)
	local neurons = {}
	if not nonInput then
		for i=1, inputs do
			neurons[i] = true
		end
	end
	for o = 1, outputs do
		neurons[maxNodes+o] = true
	end
	for i=1,#genes do
		if (not nonInput) or genes[i].into > inputs then
			neurons[genes[i].into] = true
		end
		if (not nonInput) or genes[i].out > inputs then
			neurons[genes[i].out] = true
		end
	end

	local count = 0
	for _,_ in pairs(neurons) do
		count = count + 1
	end
	local n = math.random(1, count)
	
	for k,v in pairs(neurons) do
		n = n-1
		if n == 0 then
			return k
		end
	end
	
	return 0
end

function new()
	local genome = {}
	genome.genes = {}
	genome.fitness = 0
	genome.adjustedFitness = 0
	genome.network = {}
	genome.maxneuron = 0
	genome.globalRank = 0
	genome.mutationRates = {}
	genome.mutationRates["connections"] = mutateConnectionsChance
	genome.mutationRates["link"] = linkMutationChance
	genome.mutationRates["bias"] = biasMutationChance
	genome.mutationRates["node"] = nodeMutationChance
	genome.mutationRates["enable"] = enableMutationChance
	genome.mutationRates["disable"] = disableMutationChance
	genome.mutationRates["step"] = stepSize
	
	return genome
end

function mutate(genome)
	for mutation,rate in pairs(genome.mutationRates) do
		if math.random(1,2) == 1 then
			genome.mutationRates[mutation] = 0.95*rate
		else
			genome.mutationRates[mutation] = 1.05263*rate
		end
	end

	if math.random() < genome.mutationRates["connections"] then
		genome = pointMutate(genome)
	end
	
	local p = genome.mutationRates["link"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, false)
		end
		p = p - 1
	end

	p = genome.mutationRates["bias"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, true)
		end
		p = p - 1
	end
	
	p = genome.mutationRates["node"]
	while p > 0 do
		if math.random() < p then
			nodeMutate(genome)
		end
		p = p - 1
	end
	
	p = genome.mutationRates["enable"]
	while p > 0 do
		if math.random() < p then
			enableDisableMutate(genome, true)
		end
		p = p - 1
	end

	p = genome.mutationRates["disable"]
	while p > 0 do
		if math.random() < p then
			enableDisableMutate(genome, false)
		end
		p = p - 1
	end
end

function newGenome()
	local genome = {}
	genome.genes = {}
	genome.fitness = 0
	genome.adjustedFitness = 0
	genome.network = {}
	genome.maxneuron = 0
	genome.globalRank = 0
	genome.mutationRates = {}
	genome.mutationRates["connections"] = mutateConnectionsChance
	genome.mutationRates["link"] = linkMutationChance
	genome.mutationRates["bias"] = biasMutationChance
	genome.mutationRates["node"] = nodeMutationChance
	genome.mutationRates["enable"] = enableMutationChance
	genome.mutationRates["disable"] = disableMutationChance
	genome.mutationRates["step"] = stepSize
	
	return genome
end

function mutate(genome)
	for mutation, rate in pairs(genome.mutationRates) do
		if math.random(1, 2) == 1 then
			genome.mutationRates[mutation] = 0.95 * rate
		else
			genome.mutationRates[mutation] = 1.05263 * rate
		end
	end

	if math.random() < genome.mutationRates["connections"] then
		pointMutate(genome)
	end
	
	local p = genome.mutationRates["link"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, false)
		end
		p = p - 1
	end

	p = genome.mutationRates["bias"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, true)
		end
		p = p - 1
	end
	
	p = genome.mutationRates["node"]
	while p > 0 do
		if math.random() < p then
			nodeMutate(genome)
		end
		p = p - 1
	end
	
	p = genome.mutationRates["enable"]
	while p > 0 do
		if math.random() < p then
			enableDisableMutate(genome, true)
		end
		p = p - 1
	end

	p = genome.mutationRates["disable"]
	while p > 0 do
		if math.random() < p then
			enableDisableMutate(genome, false)
		end
		p = p - 1
	end
end

function enableDisableMutate(genome, enable)
	local candidates = {}
	for _,gene in pairs(genome.genes) do
		if gene.enabled == not enable then
			table.insert(candidates, gene)
		end
	end
	
	if #candidates == 0 then
		return
	end
	
	local gene = candidates[math.random(1,#candidates)]
	gene.enabled = not gene.enabled
end

function nodeMutate(genome)
	if #genome.genes == 0 then
		return
	end

	genome.maxneuron = genome.maxneuron + 1

	local gene = genome.genes[math.random(1,#genome.genes)]
	if not gene.enabled then
		return
	end
	gene.enabled = false
	
	local gene1 = copyGene(gene)
	gene1.out = genome.maxneuron
	gene1.weight = 1.0
	gene1.innovation = newInnovation()
	gene1.enabled = true
	table.insert(genome.genes, gene1)
	
	local gene2 = copyGene(gene)
	gene2.into = genome.maxneuron
	gene2.innovation = newInnovation()
	gene2.enabled = true
	table.insert(genome.genes, gene2)
end

function basicGenome()
	local genome = newGenome()
	local innovation = 1

	genome.maxneuron = inputs
	mutate(genome)
	
	return genome
end

function disjoint(genes1, genes2)
	local i1 = {}
	for i = 1,#genes1 do
		local gene = genes1[i]
		i1[gene.innovation] = true
	end

	local i2 = {}
	for i = 1,#genes2 do
		local gene = genes2[i]
		i2[gene.innovation] = true
	end
	
	local disjointGenes = 0
	for i = 1,#genes1 do
		local gene = genes1[i]
		if not i2[gene.innovation] then
			disjointGenes = disjointGenes+1
		end
	end
	
	for i = 1,#genes2 do
		local gene = genes2[i]
		if not i1[gene.innovation] then
			disjointGenes = disjointGenes+1
		end
	end
	
	local n = math.max(#genes1, #genes2)
	
	return disjointGenes / n
end

function sameSpecies(genome1, genome2)
	local dd = deltaDisjoint*disjoint(genome1.genes, genome2.genes)
	local dw = deltaWeights*weights(genome1.genes, genome2.genes) 
	return dd + dw < deltaThreshold
end

function copyGene(gene)
	local gene2 = newGene()
	gene2.into = gene.into
	gene2.out = gene.out
	gene2.weight = gene.weight
	gene2.enabled = gene.enabled
	gene2.innovation = gene.innovation
	
	return gene2
end


function weights(genes1, genes2)
	local i2 = {}
	for i = 1,#genes2 do
		local gene = genes2[i]
		i2[gene.innovation] = gene
	end

	local sum = 0
	local coincident = 0
	for i = 1,#genes1 do
		local gene = genes1[i]
		if i2[gene.innovation] ~= nil then
			local gene2 = i2[gene.innovation]
			sum = sum + math.abs(gene.weight - gene2.weight)
			coincident = coincident + 1
		end
	end
	
	return sum / coincident
end

function newSpecies()
	local species = {}
	species.topFitness = 0
	species.staleness = 0
	species.genomes = {}
	species.averageFitness = 0
	
	return species
end

function addToSpecies(child)
	local foundSpecies = false
	for s=1,#pool.species do
		local species = pool.species[s]
		if not foundSpecies and sameSpecies(child, species.genomes[1]) then
			table.insert(species.genomes, child)
			foundSpecies = true
		end
	end
	
	if not foundSpecies then
		local childSpecies = newSpecies()
		table.insert(childSpecies.genomes, child)
		table.insert(pool.species, childSpecies)
	end
end

function generateNetwork(genome)
	local network = {}
	network.neurons = {}
	
	for i = 1, inputs do
		network.neurons[i] = newNeuron()
	end
	
	for o = 1, outputs do
		network.neurons[maxNodes+o] = newNeuron()
	end
	
	table.sort(genome.genes, function (a, b)
		return (a.out < b.out)
    end)

	for i = 1, #genome.genes do
		local gene = genome.genes[i]
		if gene.enabled then
			if network.neurons[gene.out] == nil then
				network.neurons[gene.out] = newNeuron()
			end
			local neuron = network.neurons[gene.out]
			table.insert(neuron.incoming, gene)
			if network.neurons[gene.into] == nil then
				network.neurons[gene.into] = newNeuron()
			end
		end
	end
	
    genome.network = network
    return network
end

function evaluateNetwork(network, neuralInputs)
	table.insert(neuralInputs, 1)
    if #neuralInputs ~= inputs then
        print(#neuralInputs)
        print(inputs)
		console.writeline("Incorrect number of neural network inputs.")
		return {}
	end
	
	for i = 1, inputs do
		network.neurons[i].value = neuralInputs[i]
	end
	
	for _, neuron in pairs(network.neurons) do
		local sum = 0
		for j = 1, #neuron.incoming do
			local incoming = neuron.incoming[j]
			local other = network.neurons[incoming.into]
			sum = sum + incoming.weight * other.value
		end
		
		if #neuron.incoming > 0 then
			neuron.value = sigmoid(sum)
		end
	end
	
    local neuralOutputs = {}

	for o = 1, outputs do
		local button = "P1 " .. buttonNames[o]
		if network.neurons[maxNodes+o].value > 0 then
			neuralOutputs[button] = true
		else
			neuralOutputs[button] = false
		end
	end
	
	return neuralOutputs
end

function sigmoid(x)
	return 2/(1+math.exp(-4.9*x))-1
end

function fitnessAlreadyMeasured()
	local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]
	
	return genome.fitness ~= 0
end

function nextGenome()
	pool.currentGenome = pool.currentGenome + 1
	if pool.currentGenome > #pool.species[pool.currentSpecies].genomes then
		pool.currentGenome = 1
		pool.currentSpecies = pool.currentSpecies+1
		if pool.currentSpecies > #pool.species then
			newGeneration()
			pool.currentSpecies = 1
		end
	end
end

function rankGlobally()
	local global = {}
	for s = 1,#pool.species do
		local species = pool.species[s]
		for g = 1,#species.genomes do
			table.insert(global, species.genomes[g])
		end
	end
	table.sort(global, function (a,b)
		return (a.fitness < b.fitness)
	end)
	
	for g=1,#global do
		global[g].globalRank = g
	end
end

function cullSpecies(cutToOne)
	for s = 1, #pool.species do
		local species = pool.species[s]
		
		table.sort(species.genomes, function (a,b)
			return (a.fitness > b.fitness)
		end)
		
		local remaining = math.ceil(#species.genomes/2)
		if cutToOne then
			remaining = 1
		end
		while #species.genomes > remaining do
			table.remove(species.genomes)
		end
	end
end

function removeStaleSpecies()
	local survived = {}

	for s = 1,#pool.species do
		local species = pool.species[s]
		
		table.sort(species.genomes, function (a,b)
			return (a.fitness > b.fitness)
		end)
		
		if species.genomes[1].fitness > species.topFitness then
			species.topFitness = species.genomes[1].fitness
			species.staleness = 0
		else
			species.staleness = species.staleness + 1
		end
		if species.staleness < staleSpecies or species.topFitness >= pool.maxFitness then
			table.insert(survived, species)
		end
	end

	pool.species = survived
end

function calculateAverageFitness(species)
	local total = 0
	
	for g=1,#species.genomes do
		local genome = species.genomes[g]
		total = total + genome.globalRank
	end
	
	species.averageFitness = total / #species.genomes
end

function copyGenome(genome)
	local genome2 = newGenome()
	for g=1,#genome.genes do
		table.insert(genome2.genes, copyGene(genome.genes[g]))
	end
	genome2.maxneuron = genome.maxneuron
	genome2.mutationRates["connections"] = genome.mutationRates["connections"]
	genome2.mutationRates["link"] = genome.mutationRates["link"]
	genome2.mutationRates["bias"] = genome.mutationRates["bias"]
	genome2.mutationRates["node"] = genome.mutationRates["node"]
	genome2.mutationRates["enable"] = genome.mutationRates["enable"]
	genome2.mutationRates["disable"] = genome.mutationRates["disable"]
	
	return genome2
end

function breedChild(species)
	local child = {}
	if math.random() < crossoverChance then
		g1 = species.genomes[math.random(1, #species.genomes)]
		g2 = species.genomes[math.random(1, #species.genomes)]
		child = crossover(g1, g2)
	else
		g = species.genomes[math.random(1, #species.genomes)]
		child = copyGenome(g)
	end
	
	mutate(child)
	
	return child
end

function totalAverageFitness()
	local total = 0
	for s = 1,#pool.species do
		local species = pool.species[s]
		total = total + species.averageFitness
	end

	return total
end

function removeWeakSpecies()
	local survived = {}

	local sum = totalAverageFitness()
	for s = 1,#pool.species do
		local species = pool.species[s]
		breed = math.floor(species.averageFitness / sum * population)
		if breed >= 1 then
			table.insert(survived, species)
		end
	end

	pool.species = survived
end

function newGeneration()
	cullSpecies(false) -- Cull the bottom half of each species
	rankGlobally()
	removeStaleSpecies()
	rankGlobally()
	for s = 1,#pool.species do
		local species = pool.species[s]
		calculateAverageFitness(species)
	end
	removeWeakSpecies()
	local sum = totalAverageFitness()
	local children = {}
	for s = 1,#pool.species do
		local species = pool.species[s]
		breed = math.floor(species.averageFitness / sum * population) - 1
		for i=1,breed do
			table.insert(children, breedChild(species))
		end
	end
	cullSpecies(true) -- Cull all but the top member of each species
	while #children + #pool.species < population do
		local species = pool.species[math.random(1, #pool.species)]
		table.insert(children, breedChild(species))
	end
	for c=1,#children do
		local child = children[c]
		addToSpecies(child)
	end
	
	pool.generation = pool.generation + 1

end

function writeFile(filename)
    local file = io.open(filename, "w")
file:write(pool.generation .. "\n")
file:write(pool.maxFitness .. "\n")
file:write(#pool.species .. "\n")
    for n,species in pairs(pool.species) do
    file:write(species.topFitness .. "\n")
    file:write(species.staleness .. "\n")
    file:write(#species.genomes .. "\n")
    for m,genome in pairs(species.genomes) do
        file:write(genome.fitness .. "\n")
        file:write(genome.maxneuron .. "\n")
        for mutation,rate in pairs(genome.mutationRates) do
            file:write(mutation .. "\n")
            file:write(rate .. "\n")
        end
        file:write("done\n")
        
        file:write(#genome.genes .. "\n")
        for l,gene in pairs(genome.genes) do
            file:write(gene.into .. " ")
            file:write(gene.out .. " ")
            file:write(gene.weight .. " ")
            file:write(gene.innovation .. " ")
            if(gene.enabled) then
                file:write("1\n")
            else
                file:write("0\n")
            end
        end
    end
    end
    file:close()
end

-- Diff com bonus para menos do boss
function getDiffHP()
    megaHP = tonumber(memory.readbyte(0xC2E))
    --bossHP = tonumber(memory.readbyte(0x19EE))
	--return megaHP * 0.5 - bossHP * 3 + 1000
	return megaHP
end

------------------------------------------------
------------------------------------------------
------------------------------------------------
----------------PROGRAMA------------------------
------------------------------------------------
------------------------------------------------

--draw megaman
megamanPos = getPosition()
megaX = megamanPos.x/scaleX + 21
megaY = megamanPos.y/scaleY - 3

--draw boss
bossX = memory.read_u16_le(0X19C4)
bossY = memory.readbyte(0x19C8)
bossX = bossX/scaleX + 21
bossY = bossY/scaleY + 3

initializePool(megaX, megaY, bossX, bossY)

while true do
    --draw megaman
    megamanPos = getPosition()
    megaX = megamanPos.x/scaleX + 21
    megaY = megamanPos.y/scaleY - 3
    gui.drawBox(megaX, megaY, megaX + 4, megaY + 7,0x00000000,0x80FF0000)

    --draw boss
    bossX = memory.read_u16_le(0X19C4)
    bossY = memory.readbyte(0x19C8)
    bossX = bossX/scaleX + 21
    bossY = bossY/scaleY + 3
    gui.drawBox(bossX, bossY, bossX + 4, bossY + 7,0x00000000,0x80FF0000)
    gui.drawLine(megaX, megaY + 3, bossX, bossY + 3, 0x80FF0000)

    
    local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]
    
    drawMap(genome)

	if pool.currentFrame%5 == 0 then
		evaluateCurrent(megaX, megaY, bossX, bossY)
    end
    
    joypad.set(controller)

	megaHP = tonumber(memory.readbyte(0xC2E))
	if megaHP > 20 then
		survivalMost = survival
		timeout = timeoutConstant
    end
    
	survival = survival + 1
	timeout = timeout - 1

	if timeout <= 0 then
		genome.fitness = survival
		if genome.fitness > pool.maxFitness then
            pool.maxFitness = genome.fitness
            writeFile("backup." .. pool.generation .. "." .. "megaman")
		end
        
		console.writeline("Gen " .. pool.generation .. " species " .. pool.currentSpecies .. " genome " .. pool.currentGenome .. " fitness: " .. genome.fitness)
		pool.currentSpecies = 1
        pool.currentGenome = 1
        while fitnessAlreadyMeasured() do
			nextGenome()
        end
		initializeRun(megaX, megaY, bossX, bossY)
	end

	local measured = 0
	local total = 0
	for _,species in pairs(pool.species) do
		for _,genome in pairs(species.genomes) do
			total = total + 1
			if genome.fitness ~= 0 then
				measured = measured + 1
			end
		end
    end
    
    gui.drawText(0, 0, "Gen " .. pool.generation .. " species " .. pool.currentSpecies .. " genome " .. pool.currentGenome .. " (" .. math.floor(measured/total*100) .. "%)", 0xFFFFFFFF, 11)
    gui.drawText(0, 12, "Fitness: " .. survivalMost, 0xFFFFFFFF, 11)
    gui.drawText(100, 30, "Score: " .. survival, 0xFFFFFFFF, 11)
    gui.drawText(100, 12, "Max Fitness: " .. math.floor(pool.maxFitness), 0xFFFFFFFF, 11)

	pool.currentFrame = pool.currentFrame + 1

	emu.frameadvance();
end