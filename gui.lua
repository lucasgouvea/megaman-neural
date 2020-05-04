local Gui = {}

Gui.displayGenome()
    local cells = {}
    local cell = {}
    local i = 1
    for dy = -Neural.boxRadius, Neural.boxRadius do
        for dx = -Neural.boxRadius, Neural.boxRadius do
            cell = {}
            cell.x = 60 + 5 * dx
            cell.y = 45 + 5 * dy
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
    cells[Neural.inputs] = biasCell

    BoxRadius = 6
    gui.drawBox(60-BoxRadius*5-3,45-BoxRadius*5-3,60+BoxRadius*5+2,45+BoxRadius*5+2,0xFF000000, 0x80808080)
    for n, cell in pairs(cells) do
        if n > Neural.inputs or cell.value ~= 0 then
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

    return
end 


Gui.display(pool, measured, total, rightmost, timeout, timeoutBonus)
    gui.drawText(0, 0, "Gen " .. pool.generation .. " species " .. pool.currentSpecies .. " genome " .. pool.currentGenome .. " (" .. math.floor(measured/total*100) .. "%)", 0xFF000000, 11)
    gui.drawText(0, 12, "Fitness: " .. math.floor(rightmost - (pool.currentFrame) / 2 - (timeout + timeoutBonus)*2/3), 0xFF000000, 11)
    gui.drawText(100, 12, "Max Fitness: " .. math.floor(pool.maxFitness), 0xFF000000, 11)
    return

return Gui