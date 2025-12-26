#!/bin/bash

# Auto-shutdown script for EC2 after STAR indexing completes
LOG_FILE="$HOME/rna-seq-project/auto_shutdown.log"
WORK_DIR="$HOME/rna-seq-project/work"
CHECK_INTERVAL=60  # Check every 60 seconds

echo "==================================================" | tee -a "$LOG_FILE"
echo "Auto-shutdown monitor started: $(date)" | tee -a "$LOG_FILE"
echo "==================================================" | tee -a "$LOG_FILE"

# Function to check if STAR index is complete
check_star_complete() {
    # Find the most recent STAR index directory
    STAR_INDEX_DIR=$(find "$WORK_DIR" -type d -name "star_index_dir" 2>/dev/null | head -1)
    
    if [ -z "$STAR_INDEX_DIR" ]; then
        echo "No STAR index directory found yet..." | tee -a "$LOG_FILE"
        return 1
    fi
    
    # Check for critical files that indicate completion
    if [ -f "$STAR_INDEX_DIR/Genome" ] && \
       [ -f "$STAR_INDEX_DIR/SA" ] && \
       [ -f "$STAR_INDEX_DIR/genomeParameters.txt" ]; then
        echo "✓ STAR index appears complete!" | tee -a "$LOG_FILE"
        return 0
    else
        echo "STAR indexing still in progress..." | tee -a "$LOG_FILE"
        return 1
    fi
}

# Function to verify no STAR processes are running
check_star_process() {
    if pgrep -x "STAR" > /dev/null; then
        echo "STAR process still running..." | tee -a "$LOG_FILE"
        return 1
    else
        echo "✓ No STAR processes detected" | tee -a "$LOG_FILE"
        return 0
    fi
}

# Main monitoring loop
while true; do
    echo "----------------------------------------" | tee -a "$LOG_FILE"
    echo "Check at: $(date)" | tee -a "$LOG_FILE"
    
    # Check if STAR process is running
    if pgrep -x "STAR" > /dev/null; then
        STAR_PID=$(pgrep -x "STAR")
        echo "STAR is running (PID: $STAR_PID)" | tee -a "$LOG_FILE"
        sleep "$CHECK_INTERVAL"
        continue
    fi
    
    # STAR process not running - check if index is complete
    echo "STAR process not detected, checking for completion..." | tee -a "$LOG_FILE"
    
    if check_star_complete; then
        # Double-check after 2 minutes to ensure stability
        echo "Waiting 2 minutes to confirm completion..." | tee -a "$LOG_FILE"
        sleep 120
        
        if check_star_complete && check_star_process; then
            echo "==================================================" | tee -a "$LOG_FILE"
            echo "✓ STAR indexing confirmed complete!" | tee -a "$LOG_FILE"
            echo "✓ Initiating shutdown in 60 seconds..." | tee -a "$LOG_FILE"
            echo "==================================================" | tee -a "$LOG_FILE"
            
            # Final countdown
            for i in {60..1}; do
                echo "Shutting down in $i seconds... (Ctrl+C to cancel)" | tee -a "$LOG_FILE"
                sleep 1
            done
            
            # Stop the instance
            echo "Stopping EC2 instance now: $(date)" | tee -a "$LOG_FILE"
            sudo shutdown -h now
            exit 0
        fi
    fi
    
    sleep "$CHECK_INTERVAL"
done
