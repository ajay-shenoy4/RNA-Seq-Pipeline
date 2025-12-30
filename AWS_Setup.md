# AWS Setup & Cost Analysis

## Instance Configuration

**Instance Type**: r6i.2xlarge
- vCPUs: 8
- Memory: 64 GB
- Storage: 300 GB EBS (gp3)
- Region: us-east-1

## Cost Breakdown

| Component | Cost |
|-----------|------|
| Instance (r6i.2xlarge) | $0.504/hour |
| EBS Storage (300 GB) | $0.024/hour |
| **Total** | **$0.528/hour** |

## Pipeline Runtime

- **Total**: ~45 minutes
- **Cost**: ~$0.40 per run

## Steps

1. Launch instance
2. Install Docker & Nextflow
3. Download data (~15 min)
4. Run pipeline (~45 min)
5. Download results
6. **Stop instance** (important!)

## Optimization Tips

- Use spot instances (60-70% savings)
- Pre-build AMI with tools installed
- Use S3 for data storage
- Auto-shutdown when complete
